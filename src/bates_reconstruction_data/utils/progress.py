from typing import TypeVar, Iterator, Iterable, Optional
import time


_T = TypeVar("_T")


class _ProgressIterator(Iterator[_T]):

    def __init__(self, iterable: Iterable[_T], total: Optional[int] = None, freq: float = 0.1):
        self._n: Optional[int] = total
        if self._n is None and hasattr(iterable, "__len__"):
            self._n = len(iterable)
        self._counter = 0
        self._iter = iter(iterable)
        self._start_time: Optional[float] = None
        self._last_update: float = time.monotonic()
        self._freq = freq

    def __next__(self) -> _T:
        if self._counter == 0:
            self.start()
        elif (time.monotonic() - self._last_update) > self._freq:
            self.update(self._counter, self._n)
            
        try:
            next_item = next(self._iter)
        except StopIteration:
            self.finish(self._counter, self._n)
            raise
        else:
            self._counter += 1
            return next_item

    def start(self):
        self._start_time = time.monotonic()
        self.update(0, self._n)
    
    def update(self, i: int, n: Optional[int]):
        now = time.monotonic()
        self._last_update = now

        tstr = "[    ]"
        if i > 0:
            elapsed = now - self._start_time
            
            if n is None or i == n:
                t = elapsed
            else:
                t = (elapsed / i) * (n - i)
            
            minutes = int(t // 60)
            seconds = round(t % 60)
            tstr = f"[{minutes}:{seconds:02d}]"

        frac = i / n if n is not None else None
        barlen = 30
        nbar = round(frac * barlen)
        symbol = "="
        msg = f"\r[{symbol*nbar:<{barlen}}] {frac:>5.1%} {tstr:>6}"
        print(msg, end="", flush=True)

    def finish(self, i: int, n: Optional[int]):
        self.update(i, n)
        print() # Add newline at end of progress bar


def track(iterable: Iterable[_T], total: Optional[int] = None) -> _ProgressIterator[_T]:
    return _ProgressIterator(iterable, total=total)
