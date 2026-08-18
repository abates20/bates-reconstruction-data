import functools
import time
from collections import deque


def rate_limit(n: int, t: int | float):
    """
    Rate limit a function.

    Parameters
    ----------
    n : int
        The number of allowed calls within the specified amount of time.
    t : int | float
        The amount of time in seconds during which `n` calls are allowed.
    """

    def decorator(func):
        call_queue = deque(maxlen=n)

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            if len(call_queue) >= n:
                next = call_queue.popleft() + t
                now = time.monotonic()
                call_queue.append(max(next, now))
                time.sleep(max(next - now, 0))
            else:
                call_queue.append(time.monotonic())

            return func(*args, **kwargs)
    
        return wrapper
    return decorator
