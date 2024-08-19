use clap::ValueEnum;
use std::fmt::Debug;

#[derive(Debug, PartialEq, Clone, Copy, ValueEnum)]
pub enum InputFormat {
    /// A cgmlst table with allele numbers. Optimized for ChewBBACA output
    Cgmlst,
    /// A cgmlst table with SHA1 hashes of the nucleotide of the alleles
    CgmlstHash,
    /// An alignment of nucleotide sequences in FASTA format
    Fasta,
    /// An alignment of nucleotide sequences in FASTA format. Counts all differences and not just [ACTG]
    FastaAll,
}

pub type InputMatrix = Vec<(String, SupportedTypeVec)>;

#[derive(Debug, PartialEq, Clone)]
pub enum SupportedTypeVec {
    Nucleotide(Vec<Nucleotide>),
    NucleotideAll(Vec<NucleotideAll>),
    Cgmlst(Vec<ChewBBACAinteger>),
    SHA1Hash(Vec<SHA1Hash>),
}

// #[derive(Debug, PartialEq, Clone, Copy)]
// pub enum SupportedType<T> {
//     Value(T),
// }

// impl<T> SupportedType<T>
// where
//     T: FromStr + std::fmt::Debug,
// {
//     pub fn from_str(s: &str) -> Result<Self, Box<dyn std::error::Error>>
//     where
//         <T as FromStr>::Err: 'static + std::error::Error,
//     {
//         let value = T::from_str(s)?;
//         Ok(SupportedType::Value(value))
//     }
// }

#[derive(Debug, Clone, Copy)]
pub struct ChewBBACAinteger(u16);

impl std::str::FromStr for ChewBBACAinteger {
    type Err = std::num::ParseIntError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let s = s.strip_prefix("INF-").unwrap_or(s);
        let value = u16::from_str(s).unwrap_or(0);
        Ok(ChewBBACAinteger(value))
    }
}

impl PartialEq for ChewBBACAinteger {
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0 || self.0 == 0 || other.0 == 0
    }
}

// A type that can be used to represent a fixed byte array
// and that can be parsed from a string of hex digits
#[derive(Debug, Clone, Copy)]
pub struct SHA1Hash([u8; 20]);

impl std::str::FromStr for SHA1Hash {
    type Err = std::num::ParseIntError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut bytes = [0u8; 20];
        let len = s.len() / 2;
        for i in 0..len {
            bytes[i] = u8::from_str_radix(&s[i * 2..i * 2 + 2], 16).unwrap_or_default();
        }
        Ok(SHA1Hash(bytes))
    }
}

impl PartialEq for SHA1Hash {
    fn eq(&self, other: &Self) -> bool {
        if self.0 == [0; 20] || other.0 == [0; 20] {
            return true;
        }
        self.0 == other.0
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Nucleotide(u8);

impl std::str::FromStr for Nucleotide {
    type Err = std::num::ParseIntError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "A" | "a" => Ok(Nucleotide(1)),
            "C" | "c" => Ok(Nucleotide(2)),
            "G" | "g" => Ok(Nucleotide(3)),
            "T" | "t" => Ok(Nucleotide(4)),
            _ => Ok(Nucleotide(0)),
        }
    }
}

impl PartialEq for Nucleotide {
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0 || self.0 == 0 || other.0 == 0
    }
}

impl From<u8> for Nucleotide {
    fn from(value: u8) -> Self {
        // Static lookup table for nucleotide values
        static LUT: [u8; 256] = {
            let mut lut = [0; 256];
            lut[b'a' as usize] = b'a';
            lut[b'c' as usize] = b'c';
            lut[b'g' as usize] = b'g';
            lut[b't' as usize] = b't';
            lut[b'A' as usize] = b'a';
            lut[b'C' as usize] = b'c';
            lut[b'G' as usize] = b'g';
            lut[b'T' as usize] = b't';
            lut
        };

        Nucleotide(LUT[value as usize])
    }
}
#[derive(Debug, Clone, Copy)]
pub struct NucleotideAll(u8);

impl std::str::FromStr for NucleotideAll {
    type Err = &'static str;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        // get first char from str and put it in
        s.chars()
            .next()
            .map(|c| Self(c as u8))
            .ok_or("NucleotideAll::from_str: could not parse first char from string")
    }
}

// implement construction from u8
impl From<u8> for NucleotideAll {
    fn from(value: u8) -> Self {
        Self(value.to_ascii_lowercase())
    }
}

impl PartialEq for NucleotideAll {
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::str::FromStr;

    #[test]
    fn test_chewbbaca_integer() {
        let x = ChewBBACAinteger::from_str("1").unwrap();
        assert_eq!(x, ChewBBACAinteger(1));
        let x = ChewBBACAinteger::from_str("INF-1").unwrap();
        assert_eq!(x, ChewBBACAinteger(1));
        let x = ChewBBACAinteger::from_str("INF-0").unwrap();
        assert_eq!(x, ChewBBACAinteger(0));
        let x = ChewBBACAinteger::from_str("INF-").unwrap();
        assert_eq!(x, ChewBBACAinteger(0));
    }

    #[test]
    fn test_sha1_hash() {
        let x = SHA1Hash::from_str("6bc8d04609de559621859873ef301f221cf5d991").unwrap();
        assert_eq!(
            x,
            SHA1Hash([
                0x6b, 0xc8, 0xd0, 0x46, 0x09, 0xde, 0x55, 0x96, 0x21, 0x85, 0x98, 0x73, 0xef, 0x30,
                0x1f, 0x22, 0x1c, 0xf5, 0xd9, 0x91
            ])
        );
    }

    #[test]
    fn test_nucleotide() {
        let x = Nucleotide::from_str("A").unwrap();
        assert_eq!(x, Nucleotide(1));
        let x = Nucleotide::from_str("C").unwrap();
        assert_eq!(x, Nucleotide(2));
        let x = Nucleotide::from_str("G").unwrap();
        assert_eq!(x, Nucleotide(3));
        let x = Nucleotide::from_str("T").unwrap();
        assert_eq!(x, Nucleotide(4));
        let x = Nucleotide::from_str("a").unwrap();
        assert_eq!(x, Nucleotide(1));
        let x = Nucleotide::from_str("c").unwrap();
        assert_eq!(x, Nucleotide(2));
        let x = Nucleotide::from_str("g").unwrap();
        assert_eq!(x, Nucleotide(3));
        let x = Nucleotide::from_str("t").unwrap();
        assert_eq!(x, Nucleotide(4));
        let x = Nucleotide::from_str("X").unwrap();
        assert_eq!(x, Nucleotide(0));
        assert_eq!(x, Nucleotide::from(1));
        assert_eq!(x, Nucleotide::from(2));
        assert_eq!(x, Nucleotide::from(3));
        assert_eq!(x, Nucleotide::from(4));
    }

    #[test]
    fn test_nucleotide_all() {
        let x = NucleotideAll::from_str("A").unwrap();
        assert_eq!(x, NucleotideAll(65));
        let x = NucleotideAll::from_str("C").unwrap();
        assert_eq!(x, NucleotideAll(67));
        let x = NucleotideAll::from_str("G").unwrap();
        assert_eq!(x, NucleotideAll(71));
        let x = NucleotideAll::from_str("T").unwrap();
        assert_eq!(x, NucleotideAll(84));
        let x = NucleotideAll::from_str("a").unwrap();
        assert_eq!(x, NucleotideAll(97));
        let x = NucleotideAll::from_str("c").unwrap();
        assert_eq!(x, NucleotideAll(99));
        let x = NucleotideAll::from_str("g").unwrap();
        assert_eq!(x, NucleotideAll(103));
        let x = NucleotideAll::from_str("t").unwrap();
        assert_eq!(x, NucleotideAll(116));
        let x = NucleotideAll::from_str("X").unwrap();
        assert_ne!(x, NucleotideAll::from_str("a").unwrap());
        assert_ne!(x, NucleotideAll::from_str("c").unwrap());
        assert_ne!(x, NucleotideAll::from_str("t").unwrap());
        assert_ne!(x, NucleotideAll::from_str("g").unwrap());
    }
}
