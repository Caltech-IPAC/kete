"""
Utilities for running the tutorials.
Tutorials are expected to take significantly more time than standard code to run,
and as a result are only run infrequently.

This file provides a CLI tool to run the tutorials, it must be run from the base
of the docs folder.

To run all tutorials:
>>> python utils.py build *

"""
import glob
import click
import os
from docutils.core import publish_doctree


@click.group()
def cli():
    pass


@cli.command()
def ls():
    """
    List tutorials contained within the 'tutorials' folder.
    """
    files = glob.glob("tutorials/*.rst")
    click.echo("Found tutorials:")
    click.echo("----------------")
    for file in files:
        if "index.rst" in file:
            continue
        click.echo(file)


@cli.command()
@click.argument("tutorial")
def build(tutorial):
    """
    Run tutorial code contained within the RST files in the 'tutorials' folder.
    """
    files = glob.glob("tutorials/" + tutorial + "*")
    files = [f for f in files if os.path.splitext(f)[1] == ".rst"]
    files = [f for f in files if 'index.rst' not in f]
    click.echo("Will build tutorials:")
    click.echo("---------------------")
    for file in files:
        click.echo(f"\t{file}")
    click.echo("")
    cont = input("Run these files? [Y/n]: ").strip().lower()
    match cont:
        case "yes" | "y" | "":
            pass
        case "no" | "n":
            click.echo("Exiting")
            return
        case _:
            click.echo("Command not recognized, exiting.")
            return
    for file in files:
        click.echo(f"Collecting code from : {file}")
        code = collect_code(file)
        click.echo("Running Code!")
        # pylint: disable-next=exec-used
        exec(code)
        click.echo("File complete")
    click.echo("Done")


def collect_code(file):
    """
    Given a RST file scrape all of the contents for code blocks of python.

    This will join all of the codeblocks together into one text string.
    """
    with open(file, "r", encoding='UTF8') as f:
        contents = f.read()
    doc = publish_doctree(contents)
    text = []
    for block in doc.findall(
        lambda x: x.tagname == "literal_block"
        and x.attributes["classes"] == ["code", "python"]
    ):
        text.append(block.astext())
    code = "\n".join(text)
    return code


if __name__ == "__main__":
    cli()
