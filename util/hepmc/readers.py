#============================================
# Extensions of the HepMC reader classes, to
# add random-access reading. Note that this is
# only really random-access for the ROOT-based
# reader, for ASCII it must scan the file, so
# the reading might be quite slow there.
#============================================
import sys, importlib
from util.hepmc.setup import HepMCSetup
from typing import TYPE_CHECKING, Any, Optional

if TYPE_CHECKING:
    setup = HepMCSetup(verbose=True)
    setup.PrepHepMC()
    python_dir = setup.GetPythonDirectory()
    if python_dir not in sys.path:
        sys.path = [python_dir] + sys.path
    from pyHepMC3 import HepMC3 as hm

_hepmc_module: Optional[Any] = None
_hepmc_root_io_module: Optional[Any] = None

def _get_hepmc():
    global _hepmc_module
    if _hepmc_module is None:
        setup = HepMCSetup(verbose=True)
        setup.PrepHepMC()
        python_dir = setup.GetPythonDirectory()
        if python_dir not in sys.path:
            sys.path = [python_dir] + sys.path
        _hepmc_module = importlib.import_module('pyHepMC3').HepMC3
    return _hepmc_module

def _get_hepmc_root_io():
    global _hepmc_root_io_module
    if _hepmc_root_io_module is None:
        setup = HepMCSetup(verbose=True)
        setup.PrepHepMC()
        python_dir = setup.GetPythonDirectory()
        if python_dir not in sys.path:
            sys.path = [python_dir] + sys.path
        _hepmc_root_io_module = importlib.import_module('pyHepMC3').rootIO.pyHepMC3rootIO.HepMC3
    return _hepmc_root_io_module

class ReaderAscii:
    def __init__(self, filename: str):
        hm = _get_hepmc()
        self._reader = hm.ReaderAscii(filename)
        self._current_event_number = 0
        self.filename = filename

    def read_event(self, event:'hm.GenEvent', event_number: Optional[int] = None):
        """
        Read an event from the file.

        Args:
            event: a HepMC3 event into which we'll read the results.
            event_number: If specified, skip to this event number before reading.
                         If None, read the next event.

        Returns:
            The event object, or None if no more events or invalid event number.
        """
        if event_number is not None:
            if event_number < 0:
                raise ValueError("Event number must be non-negative")

            # Calculate how many events to skip (can be negative for backwards)
            events_to_skip = event_number - self._current_event_number

            if events_to_skip != 0:
                if not self._reader.skip(events_to_skip):
                    # Skip failed - might be trying to go beyond bounds
                    return None
                self._current_event_number = event_number

        # Read the event
        self._reader.read_event(event)
        if not self._reader.failed():
            self._current_event_number += 1
        return

    def get_current_event_number(self) -> int:
        """Get the current event number (0-indexed)."""
        return self._current_event_number

    def skip_to(self, event_number: int) -> bool:
        """
        Skip to a specific event number without reading it.

        Args:
            event_number: The event number to skip to

        Returns:
            True if successful, False if the event number is out of bounds
        """
        if event_number < 0:
            raise ValueError("Event number must be non-negative")

        events_to_skip = event_number - self._current_event_number
        if events_to_skip == 0:
            return True

        if self._reader.skip(events_to_skip):
            self._current_event_number = event_number
            return True
        return False

    def __getattr__(self, name):
        # Delegate to wrapped reader, with special handling for methods that affect position
        attr = getattr(self._reader, name)

        if name == 'skip':
            def wrapped_skip(n_events):
                result = attr(n_events)
                if result:
                    self._current_event_number += n_events
                return result
            return wrapped_skip
        elif name == 'close':
            def wrapped_close():
                result = attr()
                self._current_event_number = 0
                return result
            return wrapped_close

        return attr

class ReaderRootTree:
    def __init__(self, filename: str):
        hm_root_io = _get_hepmc_root_io()
        self._reader = hm_root_io.ReaderRootTree(filename)
        self._current_event_number = 0
        self.filename = filename

    def read_event(self, event: 'hm.GenEvent', event_number: Optional[int] = None):
        """
        Read an event from the file.

        Args:
            event: a HepMC3 event into which we'll read the results.
            event_number: If specified, skip to this event number before reading.
                         If None, read the next event.

        Returns:
            The event object, or None if no more events or invalid event number.
        """
        if event_number is not None:
            if event_number < 0:
                raise ValueError("Event number must be non-negative")

            # Calculate how many events to skip (can be negative for backwards)
            events_to_skip = event_number - self._current_event_number

            if events_to_skip != 0:
                if not self._reader.skip(events_to_skip):
                    # Skip failed - might be trying to go beyond bounds
                    return None
                self._current_event_number = event_number

        # Read the event
        self._reader.read_event(event)
        if not self._reader.failed():
            self._current_event_number += 1
        return

    def get_current_event_number(self) -> int:
        """Get the current event number (0-indexed)."""
        return self._current_event_number

    def skip_to(self, event_number: int) -> bool:
        """
        Skip to a specific event number without reading it.

        Args:
            event_number: The event number to skip to

        Returns:
            True if successful, False if the event number is out of bounds
        """
        if event_number < 0:
            raise ValueError("Event number must be non-negative")

        events_to_skip = event_number - self._current_event_number
        if events_to_skip == 0:
            return True

        if self._reader.skip(events_to_skip):
            self._current_event_number = event_number
            return True
        return False

    def __getattr__(self, name):
        # Delegate to wrapped reader, with special handling for methods that affect position
        attr = getattr(self._reader, name)

        if name == 'skip':
            def wrapped_skip(n_events):
                result = attr(n_events)
                if result:
                    self._current_event_number += n_events
                return result
            return wrapped_skip
        elif name == 'close':
            def wrapped_close():
                result = attr()
                self._current_event_number = 0
                return result
            return wrapped_close

        return attr



