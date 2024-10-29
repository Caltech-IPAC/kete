//! Optional feature support for reading/writing parquet table
//!

use itertools::izip;
use itertools::Itertools;
use std::fs::File;

use crate::errors::{Error, KeteResult};
use crate::state::State;

use polars::prelude::*;

/// Write a collection of states to a parquet table.
pub fn write_states_parquet(states: &[State], filename: &str) -> KeteResult<()> {
    let desigs = Series::new(
        "desig".into(),
        states
            .iter()
            .map(|state| state.desig.to_string())
            .collect_vec(),
    );
    let jd = Series::new(
        "jd".into(),
        states.iter().map(|state| state.jd).collect_vec(),
    );
    let x = Series::new(
        "x".into(),
        states.iter().map(|state| state.pos[0]).collect_vec(),
    );
    let y = Series::new(
        "y".into(),
        states.iter().map(|state| state.pos[1]).collect_vec(),
    );
    let z = Series::new(
        "z".into(),
        states.iter().map(|state| state.pos[2]).collect_vec(),
    );
    let vx = Series::new(
        "vx".into(),
        states.iter().map(|state| state.vel[0]).collect_vec(),
    );
    let vy = Series::new(
        "vy".into(),
        states.iter().map(|state| state.vel[1]).collect_vec(),
    );
    let vz = Series::new(
        "vz".into(),
        states.iter().map(|state| state.vel[2]).collect_vec(),
    );
    let frame = Series::new(
        "frame".into(),
        states
            .iter()
            .map(|state| Into::<i32>::into(state.frame))
            .collect_vec(),
    );
    let center = Series::new(
        "center".into(),
        states.iter().map(|state| state.center_id).collect_vec(),
    );
    let mut df = DataFrame::new(vec![desigs, jd, x, y, z, vx, vy, vz, frame, center])
        .expect("Failed to construct dataframe");
    let file = File::create(filename).expect("could not create file");
    let _ = ParquetWriter::new(file)
        .finish(&mut df)
        .map_err(|_| Error::IOError("Failed to write to file".into()))?;
    Ok(())
}

/// Write a collection of states to a parquet table.
pub fn read_states_parquet(filename: &str) -> KeteResult<Vec<State>> {
    let r = File::open(filename)?;
    let reader = ParquetReader::new(r);
    let dataframe = reader.finish().map_err(|_| {
        Error::IOError("Failed to read contents of file as a parquet table.".into())
    })?;

    // dataframe.as_single_chunk_par();

    let cols = [
        "desig", "jd", "x", "y", "z", "vx", "vy", "vz", "frame", "center",
    ];

    if let Ok(columns) = dataframe.columns(cols) {
        let desigs: Vec<&str> = columns[0]
            .str()
            .expect("Designations are not all strings.")
            .into_no_null_iter()
            .collect();
        let jd: Vec<f64> = columns[1]
            .f64()
            .expect("JD are not all floats.")
            .into_no_null_iter()
            .collect();
        let x: Vec<f64> = columns[2]
            .f64()
            .expect("positions are not all floats.")
            .into_no_null_iter()
            .collect();
        let y: Vec<f64> = columns[3]
            .f64()
            .expect("positions are not all floats.")
            .into_no_null_iter()
            .collect();
        let z: Vec<f64> = columns[4]
            .f64()
            .expect("positions are not all floats.")
            .into_no_null_iter()
            .collect();
        let vx: Vec<f64> = columns[5]
            .f64()
            .expect("velocities are not all floats.")
            .into_no_null_iter()
            .collect();
        let vy: Vec<f64> = columns[6]
            .f64()
            .expect("velocities are not all floats.")
            .into_no_null_iter()
            .collect();
        let vz: Vec<f64> = columns[7]
            .f64()
            .expect("velocities are not all floats.")
            .into_no_null_iter()
            .collect();
        let frame: Vec<i32> = columns[8]
            .i32()
            .expect("Frames are not formatted correctly.")
            .into_no_null_iter()
            .collect();
        let centers: Vec<i64> = columns[9]
            .i64()
            .expect("Centers are not formatted correctly.")
            .into_no_null_iter()
            .collect();
        return Ok(izip!(desigs, jd, x, y, z, vx, vy, vz, frame, centers)
            .map(|(d, j, x0, y0, z0, vx0, vy0, vz0, f, c)| {
                State::new(
                    crate::state::Desig::Name(d.to_string()),
                    j,
                    [x0, y0, z0].into(),
                    [vx0, vy0, vz0].into(),
                    f.into(),
                    c,
                )
            })
            .collect());
    };
    Err(Error::IOError(
        "Parquet file is not a states file, missing columns.".into(),
    ))
}
