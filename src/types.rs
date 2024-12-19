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
    SHA1Hash(Vec<Hash>),
}

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
// The default hash size is 20 bytes corresponding to SHA1
// It still can be used for other hash sizes but for larger hashes the later bytes will be ignored
// Smaller hashes will be padded with zeros
#[derive(Debug, Clone, Copy)]
pub struct Hash([u8; 20]);

impl std::str::FromStr for Hash {
    type Err = std::num::ParseIntError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut bytes = [0u8; 20];
        let len = std::cmp::min(s.len() / 2, 20);
        for i in 0..len {
            bytes[i] = u8::from_str_radix(&s[i * 2..i * 2 + 2], 16).unwrap_or_default();
        }
        Ok(Hash(bytes))
    }
}

impl PartialEq for Hash {
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
            "G" | "g" => Ok(Nucleotide(4)),
            "T" | "t" => Ok(Nucleotide(8)),
            _ => Ok(Nucleotide(15)),
        }
    }
}

impl PartialEq for Nucleotide {
    fn eq(&self, other: &Self) -> bool {
        // this implementation is specific to the values of the Nucleotide enum
        (self.0 & other.0) != 0
    }
}

impl From<u8> for Nucleotide {
    fn from(value: u8) -> Self {
        // Static lookup table for nucleotide values
        static LUT: [Nucleotide; 256] = {
            let mut lut = [Nucleotide(15); 256];
            lut[b'a' as usize] = Nucleotide(1);
            lut[b'c' as usize] = Nucleotide(2);
            lut[b'g' as usize] = Nucleotide(4);
            lut[b't' as usize] = Nucleotide(8);
            lut[b'A' as usize] = Nucleotide(1);
            lut[b'C' as usize] = Nucleotide(2);
            lut[b'G' as usize] = Nucleotide(4);
            lut[b'T' as usize] = Nucleotide(8);
            lut
        };

        LUT[value as usize]
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
            .map(|c| Self::from(c as u8))
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
        let x = Hash::from_str("6bc8d04609de559621859873ef301f221cf5d991").unwrap();
        let empty_hash = Hash([0; 20]);

        assert_eq!(
            x,
            Hash([
                0x6b, 0xc8, 0xd0, 0x46, 0x09, 0xde, 0x55, 0x96, 0x21, 0x85, 0x98, 0x73, 0xef, 0x30,
                0x1f, 0x22, 0x1c, 0xf5, 0xd9, 0x91
            ])
        );

        assert_eq!(empty_hash, x);

        let sha256_hash =
            Hash::from_str("6bc8d04609de559621859873ef301f221cf5d9916bc8d04609de559621859873")
                .unwrap();
        // this should equal the first 20 bytes of the hash since later bytes are ignored in longer hashes
        assert_eq!(x, sha256_hash);

        let short_hash = Hash::from_str("6bc8d0").unwrap();
        let short_hash_padded = Hash::from_str("6bc8d0000000000000000000000000000000000").unwrap();

        assert_ne!(x, short_hash);
        assert_eq!(short_hash_padded, short_hash);
    }

    #[test]
    fn test_nucleotide() {
        let cap_a = Nucleotide::from_str("A").unwrap();
        let cap_c = Nucleotide::from_str("C").unwrap();
        let cap_g = Nucleotide::from_str("G").unwrap();
        let cap_t = Nucleotide::from_str("T").unwrap();
        let a = Nucleotide::from_str("a").unwrap();
        let c = Nucleotide::from_str("c").unwrap();
        let g = Nucleotide::from_str("g").unwrap();
        let t = Nucleotide::from_str("t").unwrap();
        let x = Nucleotide::from_str("X").unwrap();

        assert_eq!(cap_a, a);
        assert_eq!(cap_c, c);
        assert_eq!(cap_g, g);
        assert_eq!(cap_t, t);
        assert_eq!(cap_a, Nucleotide::from(b'a'));
        assert_eq!(cap_c, Nucleotide::from(b'c'));
        assert_eq!(cap_g, Nucleotide::from(b'g'));
        assert_eq!(cap_t, Nucleotide::from(b't'));

        assert_ne!(cap_a, cap_c);
        assert_ne!(cap_a, cap_g);
        assert_ne!(cap_a, cap_t);

        assert_eq!(x, Nucleotide::from(b'a'));
        assert_eq!(x, Nucleotide::from(b'c'));
        assert_eq!(x, Nucleotide::from(b'g'));
        assert_eq!(x, Nucleotide::from(b't'));
    }

    #[test]
    fn test_nucleotide_all() {
        let a = NucleotideAll::from_str("A").unwrap();
        let c = NucleotideAll::from_str("C").unwrap();
        let g = NucleotideAll::from_str("G").unwrap();
        let t = NucleotideAll::from_str("T").unwrap();
        let x = NucleotideAll::from_str("X").unwrap();

        assert_eq!(a, NucleotideAll::from(b'a'));
        assert_eq!(c, NucleotideAll::from(b'c'));
        assert_eq!(g, NucleotideAll::from(b'g'));
        assert_eq!(t, NucleotideAll::from(b't'));

        assert_ne!(a, c);
        assert_ne!(a, g);
        assert_ne!(a, t);

        assert_ne!(x, NucleotideAll::from(b'a'));
        assert_ne!(x, NucleotideAll::from(b'c'));
        assert_ne!(x, NucleotideAll::from(b'g'));
        assert_ne!(x, NucleotideAll::from(b't'));
    }
}
