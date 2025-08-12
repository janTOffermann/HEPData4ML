#============================================
# Extensions of the HepMC reader classes, to
# add random-access reading. Note that this is
# only really random-access for the ROOT-based
# reader, for ASCII it must scan the file, so
# the reading might be quite slow there.
#============================================
import sys, importlib
from util.hepmc.setup import HepMCSetup, prepend_to_pythonpath
from typing import TYPE_CHECKING, Any, Optional

if TYPE_CHECKING:
    setup = HepMCSetup(verbose=True)
    # setup.PrepHepMC()
    python_dir = setup.GetPythonDirectory()
    if python_dir not in sys.path:
        sys.path = [python_dir] + sys.path
    from pyHepMC3 import HepMC3 as hm

class ReaderAscii:
    def __init__(self, filename: str):
        self.setup = HepMCSetup(verbose=False)
        python_dir = self.setup.GetPythonDirectory()
        prepend_to_pythonpath(python_dir)

        from pyHepMC3 import HepMC3 as hm

        self.filename = filename
        self._reader = hm.ReaderAscii(self.filename)
        self._current_event_number = 0

    def read_event(self, event:'hm.GenEvent', event_number: Optional[int] = None):
        """
        Read an event from the file.

        Args:
            event: a HepMC3 event into which we'll read the results.
            event_number: If specified, skip to this event number before reading.
                         If None, read the next event.
        """
        # NOTE: Unclear if the "event_number" functionality will work particularly well,
        #       I don't think the Ascii reader allows for skipping backwards. -Jan
        from pyHepMC3 import HepMC3 as hm

        if event_number is not None:
            if event_number < 0:
                raise ValueError("Event number must be non-negative")

            # Calculate how many events to skip (can be negative for backwards)
            events_to_skip = event_number - self._current_event_number

            if events_to_skip > 0:
                if not self._reader.skip(events_to_skip):
                    # Skip failed - might be trying to go beyond bounds
                    return None
                self._current_event_number = event_number
            elif events_to_skip < 0: # skipping backwards doesn't work with AsciiReader, need to do some hackery -- reset the underlying reader
                self._reader.close()
                self._reader = hm.ReaderAscii(self.filename) # NOTE: This might slow things down? In general we'll want to use ROOT formats instead anyway.
                self._current_event_number = 0

                # now we can just call this function recursively, as we should find events_to_skip >= 0 due to the reset
                return self.read_event(event,event_number)

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
        elif(events_to_skip < 0): # reset the AsciiReader, call again
            self._reader.close()
            self._reader = hm.ReaderAscii(self.filename) # NOTE: This might slow things down?
            self._current_event_number = 0
            return self.skip_to(event_number)

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
        self.setup = HepMCSetup(verbose=False)
        python_dir = self.setup.GetPythonDirectory()
        prepend_to_pythonpath(python_dir)
        import pyHepMC3.rootIO.pyHepMC3rootIO.HepMC3 as hm_root_io

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
        """

        # With custom HepMC3, we can just do this
        if(event_number is None):
            self._reader.read_event(event)
            self._current_event_number += 1
        else:
            self._reader.read_event(event,event_number) # doesn't affect _current_event_number, because it doesn't affect the state of the underlying reader!

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



