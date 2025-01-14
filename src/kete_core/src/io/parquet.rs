//! Optional feature support for reading/writing parquet table
//!

use itertools::Itertools;
use std::fs::File;

use crate::errors::{Error, KeteResult};
use crate::frames::Frame;
use crate::state::State;

use polars::prelude::*;

/// Write a collection of states to a parquet table.
pub fn write_states_parquet(states: &[State], filename: &str) -> KeteResult<()> {
    let desigs = Column::new(
        "desig".into(),
        states
            .iter()
            .map(|state| state.desig.to_string())
            .collect_vec(),
    );
    let jd = Column::new(
        "jd".into(),
        states.iter().map(|state| state.jd).collect_vec(),
    );
    let x = Column::new(
        "x".into(),
        states.iter().map(|state| state.pos[0]).collect_vec(),
    );
    let y = Column::new(
        "y".into(),
        states.iter().map(|state| state.pos[1]).collect_vec(),
    );
    let z = Column::new(
        "z".into(),
        states.iter().map(|state| state.pos[2]).collect_vec(),
    );
    let vx = Column::new(
        "vx".into(),
        states.iter().map(|state| state.vel[0]).collect_vec(),
    );
    let vy = Column::new(
        "vy".into(),
        states.iter().map(|state| state.vel[1]).collect_vec(),
    );
    let vz = Column::new(
        "vz".into(),
        states.iter().map(|state| state.vel[2]).collect_vec(),
    );
    let frame = Column::new(
        "frame".into(),
        states
            .iter()
            .map(|state| Into::<i32>::into(state.frame))
            .collect_vec(),
    );
    let center = Column::new(
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

/// Read a collection of states from a parquet table.
pub fn read_states_parquet(filename: &str) -> KeteResult<Vec<State>> {
    // this reads the parquet table, then creates iterators over the contents, making
    // states by going through the iterators one at a time.
    let r = File::open(filename)?;
    let reader = ParquetReader::new(r);
    let mut dataframe = reader.finish().map_err(|_| {
        Error::IOError("Failed to read contents of file as a parquet table.".into())
    })?;
    let dataframe = dataframe.as_single_chunk_par();

    // create all the iterators, these are all type dependant, so they get special cased
    let mut desig_iter = dataframe
        .column("desig")
        .map_err(|_| Error::IOError("File doesn't contain the correct columns".into()))?
        .str()
        .expect("Designations are not all strings.")
        .into_no_null_iter();

    let mut frame_iter = dataframe
        .column("frame")
        .map_err(|_| Error::IOError("File doesn't contain the correct columns".into()))?
        .i32()
        .expect("Frames are not all ints.")
        .into_no_null_iter();

    let mut center_iter = dataframe
        .column("center")
        .map_err(|_| Error::IOError("File doesn't contain the correct columns".into()))?
        .i64()
        .expect("Centers are not all ints.")
        .into_no_null_iter();

    // the remaining columns are all floats, so here we make a vector of iterators of
    // floats
    let mut state_iters = dataframe
        .columns(["jd", "x", "y", "z", "vx", "vy", "vz"])
        .map_err(|_| Error::IOError("File doesn't contain the correct columns".into()))?
        .iter()
        .map(|s| {
            s.f64()
                .expect("state information is not all floats.")
                .into_no_null_iter()
        })
        .collect::<Vec<_>>();

    // throw in an assert, hopefully the compiler can see that the dimension of this is
    // then 7, meaning the index lookups below are not required to be bounds checked.
    assert!(state_iters.len() == 7);

    Ok((0..dataframe.height())
        .map(|_| {
            let desig = desig_iter
                .next()
                .expect("should have as many iterations as rows");
            let frame: Frame = frame_iter.next().unwrap().into();
            let center_id = center_iter.next().unwrap();

            let jd = state_iters[0].next().unwrap();
            let x = state_iters[1].next().unwrap();
            let y = state_iters[2].next().unwrap();
            let z = state_iters[3].next().unwrap();
            let vx = state_iters[4].next().unwrap();
            let vy = state_iters[5].next().unwrap();
            let vz = state_iters[6].next().unwrap();

            State::new(
                crate::state::Desig::Name(desig.to_string()),
                jd,
                [x, y, z].into(),
                [vx, vy, vz].into(),
                frame,
                center_id,
            )
        })
        .collect())
}
