use std::{
    fs::{self, File, OpenOptions, Permissions},
    io::{self, Write},
    path::{Path, PathBuf},
    sync::atomic::{AtomicUsize, Ordering},
};

static NEXT_TEMP_OUTPUT_ID: AtomicUsize = AtomicUsize::new(0);
static NEXT_TEMP_INPUT_ID: AtomicUsize = AtomicUsize::new(0);

pub(super) struct InputSource {
    path: PathBuf,
    _temporary: Option<TemporaryInput>,
}

impl InputSource {
    pub(super) fn from_optional_path(path: Option<&Path>) -> io::Result<Self> {
        if let Some(path) = path.filter(|path| *path != Path::new("-")) {
            return Ok(Self {
                path: path.to_path_buf(),
                _temporary: None,
            });
        }

        let temporary = TemporaryInput::from_stdin()?;
        Ok(Self {
            path: temporary.path.clone(),
            _temporary: Some(temporary),
        })
    }

    pub(super) fn path(&self) -> &Path {
        &self.path
    }
}

struct TemporaryInput {
    path: PathBuf,
}

impl TemporaryInput {
    fn from_stdin() -> io::Result<Self> {
        loop {
            let id = NEXT_TEMP_INPUT_ID.fetch_add(1, Ordering::Relaxed);
            let path = std::env::temp_dir()
                .join(format!(".fasta-util-input-{}-{id}.tmp", std::process::id()));
            let mut options = OpenOptions::new();
            options.write(true).create_new(true);
            #[cfg(unix)]
            {
                use std::os::unix::fs::OpenOptionsExt;
                options.mode(0o600);
            }
            match options.open(&path) {
                Ok(mut file) => {
                    let temporary = Self { path };
                    let mut stdin = io::stdin().lock();
                    io::copy(&mut stdin, &mut file)?;
                    file.flush()?;
                    drop(file);
                    return Ok(temporary);
                }
                Err(error) if error.kind() == io::ErrorKind::AlreadyExists => continue,
                Err(error) => return Err(error),
            }
        }
    }
}

impl Drop for TemporaryInput {
    fn drop(&mut self) {
        let _ = fs::remove_file(&self.path);
    }
}

pub(super) struct TemporaryOutput {
    destination: PathBuf,
    temporary: PathBuf,
    file: Option<File>,
    permissions: Option<Permissions>,
    committed: bool,
}

impl TemporaryOutput {
    pub(super) fn create(destination: &Path) -> io::Result<Self> {
        let destination = match fs::canonicalize(destination) {
            Ok(path) => path,
            Err(error) if error.kind() == io::ErrorKind::NotFound => destination.to_path_buf(),
            Err(error) => return Err(error),
        };
        let permissions = match fs::metadata(&destination) {
            Ok(metadata) => Some(metadata.permissions()),
            Err(error) if error.kind() == io::ErrorKind::NotFound => None,
            Err(error) => return Err(error),
        };
        let parent = destination
            .parent()
            .filter(|parent| !parent.as_os_str().is_empty())
            .unwrap_or_else(|| Path::new("."));
        if destination.file_name().is_none() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "output path has no file name",
            ));
        }

        loop {
            let id = NEXT_TEMP_OUTPUT_ID.fetch_add(1, Ordering::Relaxed);
            let temporary = parent.join(format!(".fasta-util-{}-{id}.tmp", std::process::id()));
            match OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(&temporary)
            {
                Ok(file) => {
                    return Ok(Self {
                        destination,
                        temporary,
                        file: Some(file),
                        permissions,
                        committed: false,
                    });
                }
                Err(error) if error.kind() == io::ErrorKind::AlreadyExists => continue,
                Err(error) => return Err(error),
            }
        }
    }

    pub(super) fn take_file(&mut self) -> io::Result<File> {
        self.file.take().ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                "temporary output was already opened",
            )
        })
    }

    pub(super) fn commit(&mut self) -> io::Result<()> {
        if self.file.is_some() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "temporary output is still open",
            ));
        }
        if let Some(permissions) = self.permissions.take() {
            fs::set_permissions(&self.temporary, permissions)?;
        }
        fs::rename(&self.temporary, &self.destination)?;
        self.committed = true;
        Ok(())
    }
}

impl Drop for TemporaryOutput {
    fn drop(&mut self) {
        self.file.take();
        if !self.committed {
            let _ = fs::remove_file(&self.temporary);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::TemporaryOutput;
    use std::{
        fs,
        io::Write,
        path::PathBuf,
        sync::atomic::{AtomicUsize, Ordering},
    };

    static NEXT_TEMP_DIR_ID: AtomicUsize = AtomicUsize::new(0);

    struct TemporaryDirectory(PathBuf);

    impl TemporaryDirectory {
        fn new() -> Self {
            loop {
                let id = NEXT_TEMP_DIR_ID.fetch_add(1, Ordering::Relaxed);
                let path = std::env::temp_dir().join(format!(
                    "fasta-util-output-test-{}-{id}",
                    std::process::id()
                ));
                match fs::create_dir(&path) {
                    Ok(()) => return Self(path),
                    Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => continue,
                    Err(error) => panic!("failed to create temporary directory: {error}"),
                }
            }
        }

        fn output_path(&self) -> PathBuf {
            self.0.join("result.fasta")
        }
    }

    impl Drop for TemporaryDirectory {
        fn drop(&mut self) {
            let _ = fs::remove_dir_all(&self.0);
        }
    }

    fn write_temporary_output(output: &mut TemporaryOutput, contents: &[u8]) {
        let mut file = output.take_file().unwrap();
        file.write_all(contents).unwrap();
    }

    #[test]
    fn commit_creates_destination_with_written_contents() {
        let directory = TemporaryDirectory::new();
        let destination = directory.output_path();
        let mut output = TemporaryOutput::create(&destination).unwrap();
        write_temporary_output(&mut output, b">record\nACGT");

        output.commit().unwrap();

        assert_eq!(fs::read(destination).unwrap(), b">record\nACGT");
    }

    #[test]
    fn commit_replaces_existing_destination_only_after_writing() {
        let directory = TemporaryDirectory::new();
        let destination = directory.output_path();
        fs::write(&destination, b"old contents").unwrap();
        let mut output = TemporaryOutput::create(&destination).unwrap();
        write_temporary_output(&mut output, b"new contents");

        assert_eq!(fs::read(&destination).unwrap(), b"old contents");
        output.commit().unwrap();

        assert_eq!(fs::read(destination).unwrap(), b"new contents");
    }

    #[test]
    fn dropping_uncommitted_output_preserves_destination_and_removes_temporary_file() {
        let directory = TemporaryDirectory::new();
        let destination = directory.output_path();
        fs::write(&destination, b"old contents").unwrap();
        let mut output = TemporaryOutput::create(&destination).unwrap();
        let temporary = output.temporary.clone();
        write_temporary_output(&mut output, b"incomplete contents");

        drop(output);

        assert_eq!(fs::read(destination).unwrap(), b"old contents");
        assert!(!temporary.exists());
    }

    #[test]
    fn dropping_output_before_opening_its_file_removes_temporary_file() {
        let directory = TemporaryDirectory::new();
        let destination = directory.output_path();
        let output = TemporaryOutput::create(&destination).unwrap();
        let temporary = output.temporary.clone();

        drop(output);

        assert!(!temporary.exists());
        assert!(!destination.exists());
    }

    #[test]
    fn commit_rejects_output_while_it_owns_the_file() {
        let directory = TemporaryDirectory::new();
        let destination = directory.output_path();
        fs::write(&destination, b"old contents").unwrap();
        let mut output = TemporaryOutput::create(&destination).unwrap();
        let temporary = output.temporary.clone();

        let error = output.commit().unwrap_err();
        assert_eq!(error.kind(), std::io::ErrorKind::InvalidInput);
        drop(output);

        assert_eq!(fs::read(destination).unwrap(), b"old contents");
        assert!(!temporary.exists());
    }

    #[cfg(unix)]
    #[test]
    fn commit_preserves_existing_destination_permissions() {
        use std::os::unix::fs::PermissionsExt;

        let directory = TemporaryDirectory::new();
        let destination = directory.output_path();
        fs::write(&destination, b"old contents").unwrap();
        fs::set_permissions(&destination, fs::Permissions::from_mode(0o640)).unwrap();
        let mut output = TemporaryOutput::create(&destination).unwrap();
        write_temporary_output(&mut output, b"new contents");

        output.commit().unwrap();

        let mode = fs::metadata(destination).unwrap().permissions().mode() & 0o777;
        assert_eq!(mode, 0o640);
    }
}
