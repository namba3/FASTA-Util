use fasta_util::nucleic_acid::NUCLEIC_ACID_SET;
use rand::RngExt;
use std::io::{self, BufWriter, Write};

const DEFAULT_SIZE: usize = 10_000;

fn parse_size_args<I>(args: I) -> Result<usize, io::Error>
where
    I: IntoIterator<Item = String>,
{
    let mut args = args.into_iter();
    let size = args
        .next()
        .map(|value| {
            value.parse::<usize>().map_err(|error| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("invalid sequence size '{value}': {error}"),
                )
            })
        })
        .transpose()?
        .unwrap_or(DEFAULT_SIZE);

    if let Some(extra) = args.next() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("unexpected argument '{extra}'"),
        ));
    }

    Ok(size.max(1))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let size = parse_size_args(std::env::args().skip(1))?;

    let mut rng = rand::rng();
    let set = &NUCLEIC_ACID_SET[..16];

    let stdout = io::stdout();
    let mut output = BufWriter::new(stdout.lock());
    writeln!(output, ">TestData {size} random data")?;

    let mut remaining = size;
    let mut line = [0; 50];
    while remaining > 0 {
        let line_len = remaining.min(line.len());
        for base in &mut line[..line_len] {
            *base = set[rng.random_range(0..set.len())];
        }
        output.write_all(&line[..line_len])?;
        output.write_all(b"\n")?;
        remaining -= line_len;
    }

    output.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{DEFAULT_SIZE, parse_size_args};

    fn args(values: &[&str]) -> Vec<String> {
        values.iter().map(|value| (*value).to_owned()).collect()
    }

    #[test]
    fn size_defaults_when_no_argument_is_given() {
        assert_eq!(parse_size_args(args(&[])).unwrap(), DEFAULT_SIZE);
    }

    #[test]
    fn size_uses_the_single_positional_argument() {
        assert_eq!(parse_size_args(args(&["123"])).unwrap(), 123);
    }

    #[test]
    fn zero_size_keeps_the_existing_minimum_of_one() {
        assert_eq!(parse_size_args(args(&["0"])).unwrap(), 1);
    }

    #[test]
    fn size_rejects_invalid_and_extra_arguments() {
        assert!(parse_size_args(args(&["abc"])).is_err());
        assert!(parse_size_args(args(&["123", "extra"])).is_err());
    }
}
