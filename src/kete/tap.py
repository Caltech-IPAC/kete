from __future__ import annotations
import io
import os
from functools import lru_cache
import time
import logging
from xml.etree import ElementTree
import requests
import pandas as pd
import json
import gzip
import hashlib
from .cache import cache_path


__all__ = ["tap_column_info", "query_tap"]

IRSA_TAP_URL = "https://irsa.ipac.caltech.edu/TAP/async"


logger = logging.getLogger(__name__)


@lru_cache()
def tap_column_info(table_name, base_url=IRSA_TAP_URL, auth=None):
    """
    Retrieve the column data for a specified TAP table.

    This will return a dataframe containing the column properties of the target table.

    Parameters
    ----------
    table_name :
        The name of the table to query the columns of.
    base_url :
        The URL of the TAPS service to query, this defaults to IRSA.
    auth :
        An optional (username, password), this may be used to access restricted data.
    """
    return query_tap(
        f"""SELECT * FROM TAP_SCHEMA.columns WHERE table_name='{table_name}'""",
        base_url=base_url,
        auth=auth,
    )


def query_tap(
    query,
    upload_table=None,
    base_url=IRSA_TAP_URL,
    auth=None,
    timeout=None,
    verbose=False,
    cache=True,
    update_cache=False,
):
    """
    Query TAP service, optionally upload a table which will be included in the
    query data. The pandas dataframe table will be labeled as `my_table` and columns in
    the query can be used like so:

    .. testcode::
        :skipif: True

        import kete
        import pandas as pd

        # Column names cannot match the column names in the IRSA table you are querying
        # 0 has been added to the end of these column names to satisfy this constraint.

        data = pd.DataFrame([['foo', 56823.933738, 186.249070833, 22.8977],
                            ['bar', 55232.963786, 49.14175, 21.63811111]],
                            columns=['name', 'mjd0', 'ra0', 'dec0'])

        jd = kete.Time.from_mjd(56823.933738).jd

        # This time corresponds to this phase:
        phase = kete.wise.mission_phase_from_jd(jd)

        # Single source table on IRSA is then: phase.source_table

        # The columns of data available in this source table are
        column_information = kete.tap.tap_column_info(phase.source_table)

        # Note that lots of information is available in column_information

        # Now select which columns we want IRSA to return.
        # Using TAP_UPLOAD.my_table.name we can get back the column of data we sent
        columns_to_fetch = "TAP_UPLOAD.my_table.name, mjd, ra, dec"

        query = (f"select {columns_to_fetch} from {phase.source_table} where " +
                "CONTAINS(POINT('J2000',ra,dec)," +
                "         CIRCLE('J2000'," +
                "                TAP_UPLOAD.my_table.ra0," +
                "                TAP_UPLOAD.my_table.dec0," +
                "                0.01)" +
                "         )=1 " +
                " and (((mjd - mjd0) < 0.0001) " +
                " and ((mjd0 - mjd) < 0.0001))")

        result = kete.tap.query_tap(query, upload_table=data)

    This is a blocking operation using TAP Async queries. This submits the query,
    receives a response from the TAP service containing a URL, then queries that URL
    for job status. This continues until the job either completes or errors.

    By default, queries are cached, and any calls of this function with the same
    parameters will return the cached results. This may be disabled with the `cache`
    keywords, additionally the cache can be forcably updated using `update_cache`.

    Parameters
    ----------
    query :
        An SQL text query.
    upload_table :
        An optional pandas dataframe.
    base_url :
        The URL of the TAPS service to query, this defaults to IRSA.
    auth :
        An optional (username, password), this may be used to access restricted data.
    timeout :
        Timeout for web queries. This will raise an exception if the servers do not
        respond within this time.
    verbose :
        Print status responses as they are fetched from the TAP service.
    cache :
        Bool to indicate whether or not the query should be cached.
    update_cache :
        If `cache=True`, then this value can specify if the cache will be forceable
        updated. IE: previous query results are ignored and resubmitted to the TAP
        service.
    """
    query = AsyncTapQuery(
        query=query,
        upload_table=upload_table,
        base_url=base_url,
        auth=auth,
        timeout=timeout,
        verbose=verbose,
        cache=cache,
        update_cache=update_cache,
    )
    return query.query_blocking()


class AsyncTapQuery:
    """
    Async Tap Queries

    Parameters are the same as the `query_tap` function.

    This allows for jobs to be submitted without blocking.
    """

    def __init__(
        self,
        query,
        upload_table=None,
        base_url=IRSA_TAP_URL,
        auth=None,
        timeout=None,
        verbose=False,
        cache=True,
        update_cache=False,
    ):
        self.query = query
        self.upload_table = upload_table
        self.base_url = base_url
        self.auth = auth
        self.timeout = timeout
        self.verbose = verbose
        self.update_cache = update_cache
        self.cache = cache

        self.data = dict(FORMAT="CSV", QUERY=query)
        files = None
        if upload_table is not None:
            self.data["UPLOAD"] = "my_table,param:table.tbl"

            csv_output = io.StringIO()
            pd.DataFrame(upload_table).to_csv(csv_output, index=False)
            csv_output.seek(0)
            files = {"table.tbl": csv_output.read().encode()}
        self._files = files
        if files is not None:
            _hash = int(
                hashlib.md5(str((query, tuple(files.values()))).encode()).hexdigest(),
                16,
            )
        else:
            _hash = int(hashlib.md5(str(query).encode()).hexdigest(), 16)
        self._hash = str(abs(_hash))[:16]
        path = cache_path(subfolder="tap")
        path = os.path.join(path, f"{self._hash[:3]}")
        if not os.path.isdir(path):
            os.makedirs(path)
        job_path = os.path.join(path, f"{self._hash}.json.gz")
        resp_path = os.path.join(path, f"{self._hash}.csv.gz")

        if cache and update_cache:
            if os.path.exists(job_path):
                os.remove(job_path)
            if os.path.exists(resp_path):
                os.remove(resp_path)
        os.makedirs(path, exist_ok=True)

        if cache and not os.path.exists(job_path):
            with gzip.open(job_path, "wb") as f:
                f.write(json.dumps({"status": "NOT_SUBMITTED"}).encode())
        self.job_path = job_path
        self.resp_path = resp_path
        self._status_url = None

    def query_blocking(self):
        """
        Submit the query to the TAP service, and block until the results are returned.
        """
        start = time.time()
        status = self.query_status()
        if status == "NOT_SUBMITTED":
            self.submit()
            status = "QUEUED"

        delay = 0.05
        last_print = 0
        while status in ["QUEUED", "EXECUTING"]:
            cur_time = time.time()
            elapsed = cur_time - start
            time.sleep(delay)
            status = self.query_status()

            # Increase time between queries until there is 30 seconds between.
            # Then continue forever.
            if elapsed < 2:
                pass
            elif delay < 3:
                delay += 0.05
            elif delay < 30:
                delay += 1
            if self.verbose and abs(cur_time - last_print) > 3:
                logger.info(
                    f"TAP response ({elapsed:0.1f} sec elapsed): %s",
                    status,
                )
                last_print = cur_time
        if status != "COMPLETED":
            raise ValueError("Job Failed: ", status)

        return self.result()

    def submit(self):
        """
        Submit the job to the TAP service.

        If the job results already exist in the cache, this will skip submission.
        """
        if self.cache and os.path.exists(self.resp_path):
            if self.verbose:
                logger.info(
                    (
                        "TAP query has already been completed and saved ",
                        "to cache, not submitting.",
                    ),
                )
            return

        submit = requests.post(
            self.base_url,
            data=self.data,
            files=self._files,
            auth=self.auth,
            timeout=self.timeout,
        )
        submit.raise_for_status()

        tree = ElementTree.fromstring(submit.content.decode())
        element = tree.find("{*}results")

        urls = [v for k, v in element[0].attrib.items() if "href" in k]
        if len(urls) != 1:
            raise ValueError("Unexpected results: ", submit.content.decode())
        url = urls[0]

        self._result_url = url
        self._status_url = url.replace("results/result", "phase")

        if self.cache:
            with gzip.open(self.job_path, "rb") as f:
                status_file = json.loads(f.read().decode())
            status_file["status"] = "QUEUED"
            status_file["status_url"] = self._status_url
            status_file["result_url"] = self._result_url
            with gzip.open(self.job_path, "wb") as f:
                f.write(json.dumps(status_file).encode())

    def query_status(self):
        """
        Query the status from the TAP service. If the job results already exist in
        the cache, then this will return as COMPLETED without querying.

        Status results can have one of outcomes:
        NOT_SUBMITTED, QUEUED, EXECUTING, ERROR, COMPLETED
        """
        if self.cache and os.path.exists(self.resp_path):
            return "COMPLETED"
        if self._status_url is None:
            status = "NOT_SUBMITTED"
        else:
            status = requests.get(
                self._status_url, timeout=self.timeout, auth=self.auth
            )
            status.raise_for_status()
            status = status.content.decode().upper()

        if self.cache:
            with gzip.open(self.job_path, "rb") as f:
                status_file = json.loads(f.read().decode())
            with gzip.open(self.job_path, "wb") as f:
                status_file["status"] = status
                f.write(json.dumps(status_file).encode())
        return status

    def result(self):
        """
        Fetch the finished results from the TAP service, caching the results if
        requested.

        If the results already exist in the cache, then they are returned without
        submitting any queries to the TAP service.

        Returns a Pandas Dataframe of the query results.
        """
        if self.cache and os.path.exists(self.resp_path):
            with gzip.open(self.resp_path, "rb") as f:
                return pd.read_csv(f)

        if self.verbose:
            logger.info("Downloading results...")
        result = requests.get(self._result_url, timeout=self.timeout, auth=self.auth)
        result.raise_for_status()
        result = pd.read_csv(io.StringIO(result.text))
        if self.cache:
            with gzip.open(self.resp_path, "wb") as f:
                result.to_csv(f, index=False)
            if self.verbose:
                logger.info("Results saved to cache.")

        if self.verbose:
            logger.info("Download complete.")
        return result
