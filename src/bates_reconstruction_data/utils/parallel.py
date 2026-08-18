from typing import TypeVar, Callable, Optional
import concurrent.futures

from bates_reconstruction_data.utils import progress


T = TypeVar("T")


def apply(func: Callable[..., T], *args, threads: Optional[int] = None, show_progress: bool = True) -> list[T]:
    """
    Utility method to apply a function to a list of input arguments.

    Similar to the built in `map` function, but allows for parallel execution
    across multiple threads.
    """
    lens = [len(x) for x in args if hasattr(x, "__len__")]
    n = min(lens, default=None)

    if threads == 1:
        iterator = map(func, *args)
        if show_progress:
            iterator = progress.track(iterator, total=n)
        results = list(iterator)

    else:
        executor = concurrent.futures.ThreadPoolExecutor(threads)
        try:
            futures = {executor.submit(func, *argset): i for i, argset in enumerate(zip(*args))}
            iterator = concurrent.futures.as_completed(futures)
            if show_progress:
                iterator = progress.track(iterator, total=n)
            completed = list(iterator)

            # Sort the completed futures based on input sequence and pull out the results
            completed.sort(key=lambda x: futures[x])
            results = [f.result() for f in completed]
        
        finally:
            executor.shutdown(wait=False, cancel_futures=True)
    
    return results
