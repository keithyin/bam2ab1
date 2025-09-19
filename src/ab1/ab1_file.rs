//! For reading AB1 trace files. (Applied Biosystem's sequencing)
//! [BioPython docs](https://biopython.org/wiki/ABI_traces)
//!
//! Adapted directly from this [BioPython code](https://github.com/biopython/biopython/blob/master/Bio/SeqIO/AbiIO.py)
//!
//! We are unable to find the official format spec for AB1 files.

use std::{
    ffi::{CStr, CString},
    fs::File,
    io::{self, ErrorKind, Read, Seek, SeekFrom},
    path::Path,
};

use anyhow::Context;

use crate::pascal_str::PascalString;

const HEADER_SIZE: usize = 128;
#[derive(Debug, Default)]
pub struct Header {
    pub magic_number: String,
    pub file_version: u16,
    pub dir_entry: DirEntry,
}

impl Header {
    pub fn new(file_version: u16) -> Self {
        let mut root_dir_entry = DirEntry::default();
        root_dir_entry.tag_name = "tdir".to_string();
        root_dir_entry.tag_number = 1;
        root_dir_entry.element_size = 28;
        root_dir_entry.element_type_code = 1023;

        Self {
            magic_number: "ABIF".to_string(),
            file_version: file_version,
            dir_entry: root_dir_entry,
        }
    }

    pub fn from_bytes(bytes: &[u8; HEADER_SIZE]) -> anyhow::Result<Self> {
        let mut header = Self::default();
        header.magic_number = String::from_utf8(bytes[..4].to_vec())?;
        if !header.magic_number.eq("ABIF") {
            anyhow::bail!("not a valid magic number: {}", header.magic_number);
        }

        header.file_version = u16::from_be_bytes(bytes[4..6].try_into().context("try into error")?);
        header.dir_entry = DirEntry::from_bytes(&bytes[6..(6 + DIR_ENTRY_SIZE)])?;

        Ok(header)
    }

    pub fn to_bytes(&self) -> [u8; HEADER_SIZE] {
        let mut header_bytes = [0_u8; HEADER_SIZE];
        header_bytes[..4].copy_from_slice(self.magic_number.as_bytes());
        header_bytes[4..6].copy_from_slice(&self.file_version.to_be_bytes());
        header_bytes[6..(6 + DIR_ENTRY_SIZE)].copy_from_slice(&self.dir_entry.to_bytes());
        header_bytes
    }
}

const DIR_ENTRY_SIZE: usize = 28;

#[derive(Debug, Default, Clone)]
pub struct DirEntry {
    pub tag_name: String, // 4 bytes always.
    pub tag_number: u32,
    pub element_type_code: u16,
    pub element_size: u16,
    pub num_elements: u32,
    pub data_size: u32,
    pub data_offset: u32,
    pub data_handle: u32,
}

impl DirEntry {
    pub fn from_bytes(bytes: &[u8]) -> anyhow::Result<Self> {
        assert!(bytes.len() == DIR_ENTRY_SIZE);
        Ok(Self {
            tag_name: String::from_utf8(bytes[..4].to_vec())?,
            tag_number: u32::from_be_bytes(bytes[4..8].try_into().context("try into error")?),
            element_type_code: u16::from_be_bytes(
                bytes[8..10].try_into().context("try into error")?,
            ),
            element_size: u16::from_be_bytes(bytes[10..12].try_into().context("try into error")?),
            num_elements: u32::from_be_bytes(bytes[12..16].try_into().context("try into error")?),
            data_size: u32::from_be_bytes(bytes[16..20].try_into().context("try into error")?),
            data_offset: u32::from_be_bytes(bytes[20..24].try_into().context("try into error")?),
            data_handle: u32::from_be_bytes(bytes[24..28].try_into().context("try into error")?),
        })
    }

    pub fn to_bytes(&self) -> [u8; DIR_ENTRY_SIZE] {
        let mut result_bytes = [0_u8; DIR_ENTRY_SIZE];
        result_bytes[..4].copy_from_slice(self.tag_name.as_bytes());
        result_bytes[4..8].copy_from_slice(&self.tag_number.to_be_bytes());
        result_bytes[8..10].copy_from_slice(&self.element_type_code.to_be_bytes());
        result_bytes[10..12].copy_from_slice(&self.element_size.to_be_bytes());
        result_bytes[12..16].copy_from_slice(&self.num_elements.to_be_bytes());
        result_bytes[16..20].copy_from_slice(&self.data_size.to_be_bytes());
        result_bytes[20..24].copy_from_slice(&self.data_offset.to_be_bytes());
        result_bytes[24..28].copy_from_slice(&self.data_handle.to_be_bytes());

        result_bytes
    }
}

#[derive(Debug)]
struct AbiFileReader<R: Read + Seek> {
    stream: R,
    header: Header,
    tagged_datas: Vec<TaggedData>,
}

impl<R: Read + Seek> AbiFileReader<R> {
    pub fn new(mut stream: R) -> anyhow::Result<Self> {
        let mut header_bytes = [0_u8; HEADER_SIZE];

        stream
            .read_exact(&mut header_bytes)
            .context("read header error")?;

        let header = Header::from_bytes(&header_bytes).context("bytes to header error")?;
        eprintln!("{header:?}");

        Ok(Self {
            stream,
            header,
            tagged_datas: vec![],
        })
    }

    pub fn get_tagged_datas(&self) -> &[TaggedData] {
        &self.tagged_datas
    }

    pub fn extract_tagged_datas(&mut self) -> anyhow::Result<()> {
        for i in 0..self.header.dir_entry.num_elements {
            // todo: QC data_offset; coming out much too high.
            // Note: Element size should always be DIR_SIZE.
            let dir_entry_start = self.header.dir_entry.data_offset as usize
                + (i * self.header.dir_entry.element_size as u32) as usize;

            self.stream.seek(SeekFrom::Start(dir_entry_start as u64))?;
            let mut dir_buf = [0; DIR_ENTRY_SIZE];
            if self.stream.read(&mut dir_buf)? == 0 {
                break;
            };

            let dir = DirEntry::from_bytes(&dir_buf)?;

            let tag_buf = if dir.data_size <= 4 {
                let mut tag_buf = vec![0_u8; dir.data_size as usize];
                tag_buf.copy_from_slice(&dir.data_offset.to_be_bytes()[..dir.data_size as usize]);
                tag_buf
            } else {
                self.stream.seek(SeekFrom::Start(dir.data_offset as u64))?;
                let mut tag_buf = vec![0; dir.data_size as usize];
                if self.stream.read(&mut tag_buf)? == 0 {
                    break;
                };
                tag_buf
            };

            let tag_data =
                parse_tag_data(dir.element_type_code, dir.num_elements as usize, &tag_buf)?;
            self.tagged_datas.push(TaggedData {
                tag: dir.tag_name.clone(),
                tag_number: dir.tag_number,
                data: tag_data,
            });
        }

        return Ok(());
    }
}

pub struct AbiFile {
    header: Header,
    tagged_datas: Vec<TaggedData>,
}

impl AbiFile {
    pub fn new(version: u16) -> Self {
        let header = Header::new(version);
        Self {
            header,
            tagged_datas: vec![],
        }
    }

    pub fn push_tagged_data(&mut self, tagged_data: TaggedData) {
        self.header.dir_entry.element_size += 1;
        self.header.dir_entry.data_size += self.header.dir_entry.element_size as u32;

        self.tagged_datas.push(tagged_data);
    }
}

#[derive(Debug, Clone)]
pub struct TaggedData {
    pub tag: String, // 4 bytes string
    pub tag_number: u32,
    pub data: Data,
}

impl TaggedData {
    pub fn build_dir_entry(&self) -> DirEntry {
        let mut entry = DirEntry::default();
        entry.tag_name = self.tag.clone();
        entry.tag_number = self.tag_number;
        entry.element_type_code = self.data.element_code();
        entry.element_size = self.data.get_ele_size();
        entry.num_elements = self.data.get_num_ele();
        entry.data_size = self.data.get_data_size();

        if entry.data_size <= 4 {
            let mut b16_bytes = [0_u8; 4];
            let data_bytes = self.data.to_bytes();
            let copy_start = 4 - data_bytes.len();
            b16_bytes[copy_start..].copy_from_slice(&data_bytes);
            entry.data_offset = u32::from_be_bytes(b16_bytes);
        }

        entry
    }
}

#[derive(Debug, Clone)]
pub enum Data {
    U8(Vec<u8>),
    U16(Vec<u16>),
    U32(Vec<u32>),
    String(String),
    PascallString(PascalString),
    CString(CString),
}

impl Data {
    pub fn get_ele_size(&self) -> u16 {
        match self {
            Self::U8(_) => 1,
            Self::U16(_) => 2,
            Self::U32(_) => 4,
            Self::String(_) => 1,
            Self::PascallString(_) => 1,
            Self::CString(_) => 1,
        }
    }

    pub fn get_num_ele(&self) -> u32 {
        let num_ele = match self {
            Self::U8(v) => v.len(),
            Self::U16(v) => v.len(),
            Self::U32(v) => v.len(),
            Self::String(v) => v.len(),
            Self::PascallString(v) => v.bytes_len(),
            Self::CString(v) => v.as_bytes().len() + 1,
        };
        num_ele as u32
    }

    pub fn get_data_size(&self) -> u32 {
        self.get_ele_size() as u32 * self.get_num_ele()
    }

    pub fn to_bytes(&self) -> Vec<u8> {
        match self {
            Self::U8(v) => v.to_vec(),
            Self::U16(v) => v.iter().flat_map(|v| v.to_be_bytes().into_iter()).collect(),
            Self::U32(v) => v.iter().flat_map(|v| v.to_be_bytes().into_iter()).collect(),
            Self::String(v) => v.as_bytes().to_vec(),
            Self::PascallString(v) => v.bytes_raw(),
            Self::CString(v) => v.as_bytes_with_nul().to_vec(),
        }
    }

    pub fn element_code(&self) -> u16 {
        match self {
            Self::U8(_) => 1,
            Self::U16(_) => 4,
            Self::U32(_) => 5,
            Self::String(_) => 2,
            Self::PascallString(_) => 18,
            Self::CString(_) => 19,
        }
    }

    // used for print
    pub fn truncated_result(&self) -> String {
        let print_ele = self.get_num_ele().min(100) as usize;
        match self {
            Self::U8(v) => format!("{:?}", &v[..print_ele]),
            Self::U16(v) => format!("{:?}", &v[..print_ele]),
            Self::U32(v) => format!("{:?}", &v[..print_ele]),
            Self::String(v) => v.to_string(),
            Self::PascallString(v) => v.to_string().unwrap(),
            Self::CString(v) => v.to_str().unwrap().to_string(),
        }
    }
}

fn parse_tag_data(elem_code: u16, _elem_num: usize, data: &[u8]) -> io::Result<Data> {
    //     1: "b",  # byte
    //     2: "s",  # char
    //     3: "H",  # word
    //     4: "h",  # short
    //     5: "i",  # long
    //     6: "2i",  # rational, legacy unsupported
    //     7: "f",  # float
    //     8: "d",  # double
    //     10: "h2B",  # date
    //     11: "4B",  # time
    //     12: "2i2b",  # thumb
    //     13: "B",  # bool
    //     14: "2h",  # point, legacy unsupported
    //     15: "4h",  # rect, legacy unsupported
    //     16: "2i",  # vPoint, legacy unsupported
    //     17: "4i",  # vRect, legacy unsupported
    //     18: "s",  # pString, use the leading 1 byte to represent the length of the string
    //     19: "s",  # cString, null-terminated string
    //     20: "2i",  # tag, legacy unsupported

    match elem_code {
        2 => Ok(Data::String(
            std::str::from_utf8(data).unwrap_or("").to_string(),
        )),
        4 => {
            let as_u16 = data
                .chunks_exact(2)
                .map(|chunk| u16::from_be_bytes([chunk[0], chunk[1]]))
                .collect();
            Ok(Data::U16(as_u16))
        }
        5 => {
            let as_u32 = data
                .chunks_exact(4)
                .map(|chunk| u32::from_be_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]))
                .collect();
            Ok(Data::U32(as_u32))
        }

        18 => Ok(Data::PascallString(PascalString::new(data.to_vec()))),
        19 => Ok(Data::CString(unsafe {
            CStr::from_ptr(data.as_ptr() as *const i8).to_owned()
        })),

        _ => {
            // todo: Handle appropriately.
            Err(io::Error::new(
                ErrorKind::InvalidData,
                format!("Invalid data type in AB1 file: {elem_code}"),
            ))
        }
    }
}

/// Read a file in the GenBank format.
/// [Rust docs ref of fields](https://docs.rs/gb-io/latest/gb_io/seq/struct.Seq.html)
pub fn import_ab1(path: &Path) -> anyhow::Result<Vec<TaggedData>> {
    let file = File::open(path)?;
    let mut iterator = AbiFileReader::new(file)?;
    iterator.extract_tagged_datas()?;
    Ok(iterator.get_tagged_datas().to_vec())
}

// PLOC1: short[]
pub fn build_ploc1(data: Vec<u16>) -> TaggedData {
    TaggedData {
        tag: format!("PLOC"),
        tag_number: 1,
        data: Data::U16(data),
    }
}
// FWO_1: char[] base order sequencing analysis filter wheel order. fixed for 3500 GATC

pub fn build_fwo_1() -> TaggedData {
    TaggedData {
        tag: format!("FWO_"),
        tag_number: 1,
        data: Data::String(format!("GATC")),
    }
}
// LANE:  short[] 2bytes

pub fn build_lane() -> TaggedData {
    TaggedData {
        tag: format!("LANE"),
        tag_number: 1,
        data: Data::U16(vec![42]),
    }
}

// S/N%1: short[], signal strength for each dye

pub fn build_snr() -> TaggedData {
    TaggedData {
        tag: format!("S/N%"),
        tag_number: 1,
        data: Data::U16(vec![53, 75, 79, 48]),
    }
}

// SMPL1: pstring, sample name

pub fn build_smpl_1(name: String) -> TaggedData {
    TaggedData {
        tag: format!("SMPL"),
        tag_number: 1,
        data: Data::PascallString(
            name.as_str()
                .try_into()
                .expect("string -> pascal string error"),
        ),
    }
}

// PBAS1: char[] , array of sequence characters edited by user. PBAS2: by basecaller

pub fn build_pbas_1(bases: String) -> TaggedData {
    TaggedData {
        tag: format!("PBAS"),
        tag_number: 1,
        data: Data::String(bases),
    }
}

// CMNT1: pString, comments about sample
pub fn build_cmnt_1(comment: String) -> TaggedData {
    TaggedData {
        tag: format!("CMNT"),
        tag_number: 1,
        data: Data::PascallString(
            comment
                .as_str()
                .try_into()
                .expect("string -> pascal string error"),
        ),
    }
}
// PDMF1: pstring. sequencing analysis mobility file name chosen in collection
pub fn build_pdmf_1(fname: String) -> TaggedData {
    TaggedData {
        tag: format!("PDMF"),
        tag_number: 1,
        data: Data::PascallString(
            fname
                .as_str()
                .try_into()
                .expect("string -> pascal string error"),
        ),
    }
}

// PDMF2: pstring. ..
pub fn build_pdmf_2(fname: String) -> TaggedData {
    TaggedData {
        tag: format!("PDMF"),
        tag_number: 2,
        data: Data::PascallString(
            fname
                .as_str()
                .try_into()
                .expect("string -> pascal string error"),
        ),
    }
}
// DATA9: short array holding holding analyzed color data
pub fn build_data(tag_number: u32, data: Vec<u16>) -> TaggedData {
    TaggedData {
        tag: format!("DATA"),
        tag_number: tag_number,
        data: Data::U16(data),
    }
}
// DATA10:
// DATA11:
// DATA12:

#[cfg(test)]
mod test {
    use std::path;

    #[test]
    fn test_read_ab1() {
        let fpath = "/data-slow/kangwei-deliver/kangwei-deliver/S22509070002-Epi5A-1.ab1";
        let p = path::Path::new(fpath);
        let records = super::import_ab1(p).expect("import ab1 error");
        records.iter().for_each(|tagged_data| {
            println!(
                "{}{}:{}->{}",
                tagged_data.tag,
                tagged_data.tag_number,
                tagged_data.data.get_num_ele(),
                tagged_data.data.truncated_result()
            )
        });
        // println!("{records:?}");
    }
}
