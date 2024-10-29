//! Converting to and from bytes
use crate::errors::{Error, KeteResult};
use std::io::Read;

/// Read the exact number of specified bytes from the file.
pub(crate) fn read_bytes_exact<T: Read>(buffer: T, n_bytes: usize) -> KeteResult<Box<[u8]>> {
    let mut bytes = Vec::with_capacity(n_bytes);
    let n_read = buffer.take(n_bytes as u64).read_to_end(&mut bytes)?;
    if n_read != n_bytes {
        Err(Error::IOError("Unexpected end of file.".into()))?;
    };

    Ok(bytes.into())
}

/// Change a collection of bytes into a f64.
pub(crate) fn bytes_to_f64(bytes: &[u8], little_endian: bool) -> KeteResult<f64> {
    let bytes: [u8; 8] = bytes
        .try_into()
        .map_err(|_| Error::IOError("File is not correctly formatted".into()))?;
    if little_endian {
        Ok(f64::from_le_bytes(bytes))
    } else {
        Ok(f64::from_be_bytes(bytes))
    }
}

/// Change a collection of bytes into a vector of f64s.
pub(crate) fn bytes_to_f64_vec(bytes: &[u8], little_endian: bool) -> KeteResult<Box<[f64]>> {
    let byte_len = bytes.len();
    if byte_len % 8 != 0 {
        Err(Error::IOError("File is not correctly formatted".into()))?;
    }
    let res: Box<[f64]> = (0..byte_len / 8)
        .map(|idx| bytes_to_f64(&bytes[8 * idx..(8 + 8 * idx)], little_endian))
        .collect::<KeteResult<_>>()?;
    Ok(res)
}

/// Change a collection of bytes into a vector of i32s.
pub(crate) fn bytes_to_i32_vec(bytes: &[u8], little_endian: bool) -> KeteResult<Box<[i32]>> {
    let byte_len = bytes.len();
    if byte_len % 4 != 0 {
        Err(Error::IOError("File is not correctly formatted".into()))?;
    }
    let res: Box<[i32]> = (0..byte_len / 4)
        .map(|idx| bytes_to_i32(&bytes[4 * idx..(4 + 4 * idx)], little_endian))
        .collect::<KeteResult<_>>()?;
    Ok(res)
}

/// Change a collection of bytes into a i32.
pub(crate) fn bytes_to_i32(bytes: &[u8], little_endian: bool) -> KeteResult<i32> {
    let bytes: [u8; 4] = bytes
        .try_into()
        .map_err(|_| Error::IOError("File is not correctly formatted".into()))?;
    if little_endian {
        Ok(i32::from_le_bytes(bytes))
    } else {
        Ok(i32::from_be_bytes(bytes))
    }
}

/// Change a collection of bytes into a String.
pub(crate) fn bytes_to_string(bytes: &[u8]) -> String {
    let mut bytes = bytes.to_vec();
    bytes.iter_mut().for_each(|x| {
        if x == &0x00 {
            *x = 0x0a;
        }
    });
    String::from_utf8_lossy(&bytes).to_string()
}

/// Read a multiple contiguous f64s from the file.
pub(crate) fn read_f64_vec<T: Read>(
    buffer: T,
    n_floats: usize,
    little_endian: bool,
) -> KeteResult<Box<[f64]>> {
    let bytes = read_bytes_exact(buffer, 8 * n_floats)?;
    bytes_to_f64_vec(&bytes, little_endian)
}

/// Read a string of the specified length from the file.
/// 0x00 are replaced with new lines, and new lines are stripped from the end of the
/// string.
pub(crate) fn read_str<T: Read>(buffer: T, length: usize) -> KeteResult<String> {
    let bytes = read_bytes_exact(buffer, length)?;
    Ok(bytes_to_string(&bytes))
}
