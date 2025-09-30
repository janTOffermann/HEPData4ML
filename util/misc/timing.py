import time, functools, datetime
from collections import defaultdict
from contextlib import contextmanager
from threading import local
import numpy as np

class BasicTimer:
    def __init__(self):
        self.dict = {}

    def start_main(self):
        self.start_time = time.time()

    def end_main(self):
        self.end_time = time.time()

    # Function to help with timestamps (for our very basic profiling)
    def start_timestamp(self, key):
        if(key) not in self.dict.keys():
            self.dict[key] = {'start':[],'end':[]}
        self.dict[key]['start'].append(time.time())

    def end_timestamp(self, key):
        if(key) not in self.dict.keys():
            self.dict[key] = {'start':[],'end':[]}
        self.dict[key]['end'].append(time.time())

    def _summarize_main(self):
        elapsed_time = self.end_time - self.start_time
        elapsed_time_readable = str(datetime.timedelta(seconds=elapsed_time))
        print('Time elapsed = {:.1f} seconds.'.format(elapsed_time))
        print('({})'.format(elapsed_time_readable))

    def summarize_time(self):
        print('\n#############################')
        self._summarize_main()
        print('Breakdown by step:')
        for key in self.dict.keys():
            elapsed_time = np.sum(np.array(self.dict[key]['end']) - np.array(self.dict[key]['start']))
            elapsed_time_readable = str(datetime.timedelta(seconds=elapsed_time))
            print('\tTime on {} step: {:.1f} seconds\t({})'.format(key, elapsed_time,elapsed_time_readable))
        print('\n#############################')
        return


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
        width = int(np.max([len(x) for x in self.times.keys()]) + 3)
        print(f"{'Method':<{width}} {'Total Time (s)':<15} {'Avg Time (s)':<15} {'Calls':<10}")
        print("-" * (40 + width))
        for method, total_time in sorted(self.total_times.items(), key=lambda x: x[1], reverse=True):
            avg_time = total_time / len(self.times[method])
            calls = len(self.times[method])
            print(f"{method:<{width}} {total_time:<15.4f} {avg_time:<15.6f} {calls:<10}")

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