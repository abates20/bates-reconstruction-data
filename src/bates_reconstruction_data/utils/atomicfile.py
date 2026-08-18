from typing import Union, Literal, Optional
from tempfile import NamedTemporaryFile
import os


class AtomicFileWrapper:
    """
    A file wrapper to handle writing data to a temporary file and then moving the
    file to a specified path once all data has been written.

    For example, the intended use is for avoiding files that are partially
    downloaded by first writing the downloaded data to a temporary location and
    then moving the temporary file to the permanent location once all the data
    has been downloaded.
    """

    def __init__(self, path: Union[str, bytes, os.PathLike], mode: Literal["t", "b"] = "t"):
        self.path = str(path)
        self.mode: Literal["t", "b"] = mode
        self._file = None
        self._tmppath = None

    def __enter__(self):
        dirname, basename = os.path.split(self.path)
        self._file = NamedTemporaryFile(
            f"x+{self.mode}",
            suffix=".tmp",
            prefix="~" + basename + "-",
            dir=dirname,
            delete=False,
        )
        self._tmppath = self._file.name
        return self

    def __exit__(self, exc_type, exc_value, exc_tb):
        self.close()
        if exc_type is None:  # No error
            os.replace(self._tmppath, self.path)
        self.cleanup()
        return False

    def write(self, data: Union[str, bytes]):
        """
        Write data to the file.
        """
        self._file.write(data)

    def close(self):
        """
        Close the temporary file, flushing all data in the buffer.
        """
        if self._file is not None and not self._file.closed:
            self._file.flush()
            os.fsync(self._file.fileno())
            self._file.close()

    def cleanup(self):
        """
        Delete the temporary file if it still exists.
        """
        if self._tmppath is not None and os.path.exists(self._tmppath):
            os.remove(self._tmppath)

    def seek(self, target: int, whence: int = 0) -> int:
        return self._file.seek(target, whence)
    
    def tell(self) -> int:
        return self._file.tell()
    
    def flush(self) -> None:
        return self._file.flush()
    
    def read(self, size: Optional[int] = -1):
        return self._file.read(size)
    
    def readable(self) -> bool:
        return self._file.readable()
    
    def writable(self) -> bool:
        return self._file.writable()
    
    def seekable(self) -> bool:
        return self._file.seekable()
