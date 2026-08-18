from typing import Callable, Union, Optional, Any
import os
from tempfile import NamedTemporaryFile
import time

import requests


CallbackT = Callable[[int, int, int], Any]


def download(
    url: str,
    path: Union[str, bytes, os.PathLike],
    callback: Optional[CallbackT] = None,
    method: str = "GET",
    **kwargs
):
    """
    Download the contents of a url and save to the specified path.
    """

    # Download contents (first to temporary path and then gets renamed when exiting `with` block)
    with requests.request(method, url, stream=True, **kwargs) as response, DownloadFileObj(path) as file:
        response: requests.Response
        response.raise_for_status()
        
        total_size = int(response.headers.get("content-length", 0))
        block_size = 64 * 1024
        count = 0

        for chunk in response.iter_content(chunk_size=block_size):
            count += 1
            file.write(chunk)

            if callback is not None:
                callback(count, block_size, total_size)

    return path


class DownloadFileObj:
    """
    A file object to handle writing data to a temporary file and then moving the
    file to a specified path once all data has been written.

    For example, the intended use is for avoiding files that are partially
    downloaded by first writing the downloaded data to a temporary location and
    then moving the temporary file to the permanent location once all the data
    has been downloaded.
    """

    def __init__(self, path: Union[str, bytes, os.PathLike]):
        self.path = str(path)
        self._file = None
        self._tmppath = None

    def __enter__(self):
        dirname, basename = os.path.split(self.path)
        self._file = NamedTemporaryFile(
            "xb",
            suffix=".tmp",
            prefix=basename + "-",
            dir=dirname,
            delete=False,
        )
        self._tmppath = self._file.name
        return self

    def __exit__(self, exc_type, exc_value, exc_tb):
        self.close()
        if exc_type is None:  # Successful download
            os.replace(self._tmppath, self.path)
        self.cleanup()
        return False

    def write(self, data: bytes):
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


class DownloadProgress:
    """
    A callback object to display the progress of a download.
    """

    def __init__(self, msg: str = "Downloading...", freq: float = 0.05):
        self.msg: str = msg
        self.freq: float = freq
        self._prev: Optional[float] = None

    def __call__(self, count: int, block: int, total: int) -> None:
        now = time.monotonic()
        finished = count * block >= total
        if self._prev is None or finished or (now - self._prev) >= self.freq:
            self._prev = now
            progress = f"\r{self.msg} {min(count*block/total, 1):.1%}"
            end = "\n" if finished else ""
            print(progress, end=end, flush=True)
