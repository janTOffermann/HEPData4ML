import time
import functools
from collections import defaultdict
from contextlib import contextmanager
from threading import local

class TimingProfiler:
    def __init__(self):
        self.times = defaultdict(list)
        self.total_times = defaultdict(float)

    def profile_method(self, name=None):
        def decorator(func):
            method_name = name or f"{func.__qualname__}"

            @functools.wraps(func)
            def wrapper(*args, **kwargs):
                # Get current profiler from thread-local storage
                profiler = _get_current_profiler()
                if profiler is None:
                    return func(*args, **kwargs)

                start_time = time.perf_counter()
                result = func(*args, **kwargs)
                elapsed = time.perf_counter() - start_time

                profiler.times[method_name].append(elapsed)
                profiler.total_times[method_name] += elapsed
                return result
            return wrapper
        return decorator

    def report(self):
        print(f"{'Method':<40} {'Total Time (s)':<15} {'Avg Time (s)':<15} {'Calls':<10}")
        print("-" * 80)
        for method, total_time in sorted(self.total_times.items(), key=lambda x: x[1], reverse=True):
            avg_time = total_time / len(self.times[method])
            calls = len(self.times[method])
            print(f"{method:<40} {total_time:<15.4f} {avg_time:<15.6f} {calls:<10}")

# Thread-local storage for current profiler
_thread_local = local()

def _get_current_profiler():
    return getattr(_thread_local, 'profiler', None)

@contextmanager
def profiling_context():
    """Context manager that sets up profiling for the current thread."""
    profiler = TimingProfiler()
    old_profiler = getattr(_thread_local, 'profiler', None)
    _thread_local.profiler = profiler

    try:
        yield profiler
    finally:
        _thread_local.profiler = old_profiler

def profile_method(name=None):
    """Decorator that profiles methods when within a profiling context."""
    def decorator(func):
        method_name = name or f"{func.__qualname__}"

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            profiler = _get_current_profiler()
            if profiler is None:
                return func(*args, **kwargs)

            start_time = time.perf_counter()
            result = func(*args, **kwargs)
            elapsed = time.perf_counter() - start_time

            profiler.times[method_name].append(elapsed)
            profiler.total_times[method_name] += elapsed
            return result
        return wrapper
    return decorator

@contextmanager
def profile_block(block_name):
    """Context manager for profiling arbitrary code blocks."""
    profiler = _get_current_profiler()
    if profiler is None:
        yield
        return

    start_time = time.perf_counter()
    try:
        yield
    finally:
        elapsed = time.perf_counter() - start_time
        profiler.times[block_name].append(elapsed)
        profiler.total_times[block_name] += elapsed