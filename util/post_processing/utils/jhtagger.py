import os,glob,pathlib,threading, queue, re
from collections import deque
from typing import TYPE_CHECKING
import subprocess as sub
import numpy as np
import ROOT as rt
from util.calcs import embed_array_inplace
from util.qol_utils.progress_bar import printProgressBar, printProgressWithOutput


if TYPE_CHECKING: # Only imported during type checking -- avoids circular imports we'd otherwise get, since jets imports this file
    from util.post_processing.jets import JetFinder

class JHTaggerSetup:
    """
    This class is for setting up the Johns Hopkins top tagger (JHTagger) code.
    This is some custom C++/ROOT code, that leverages the JHTagger implementation
    in the Fastjet library.
    """
    def __init__(self,configurator,jhtagger_dir=None):
        self.SetDirectory(jhtagger_dir)
        self.SetConfigurator(configurator)
        self.executable = '/bin/bash'
        self.lib_extensions = ['so','dylib']
        self.status = False

    def SetDirectory(self,dir):
        if(dir is None):
            dir = os.path.dirname(os.path.abspath(__file__)) + '/../../root/jhtagger'
        self.dir = dir
        self.logfile = '{}/log.stdout'.format(self.dir)
        self.errfile = '{}/log.stderr'.format(self.dir)

    def SetConfigurator(self,configurator):
        self.configurator = configurator

    def FullPreparation(self):
        self.status = False
        try: self.LoadJHTagger(quiet=True)
        except:
            self.BuildJHTagger()
            self.LoadJHTagger()
        self.status = True
        return

    def BuildJHTagger(self):
        """
        Builds the JHTagger C++/ROOT library.
        This runs a shell script, where certain
        environment variables have been set by
        this Python script (via our configurator,
        which tells us where Fastjet is installed).
        """
        self.build_script = 'build.sh'
        fastjet_dir = self.configurator.GetFastjetDirectory()

        # TODO: Clean this up. Somehow put information from fastjet setup back into configurator?
        if(fastjet_dir is None):
            fastjet_dir = os.path.dirname(os.path.abspath(__file__)) + '/../../../external/fastjet'

        # TODO: This might crash if fastjet is not built yet? Could be an issue for a first run, need to test!
        fastjet_lib = glob.glob('{}/**/lib'.format(fastjet_dir),recursive=True)[0]
        fastjet_inc = glob.glob('{}/**/include'.format(fastjet_dir),recursive=True)[0]
        env = os.environ.copy()
        env['CMAKE_PREFIX_PATH'] += ':{}'.format(fastjet_lib)
        env['FASTJET_INCLUDE_DIR'] = fastjet_inc

        command = ['./{}'.format(self.build_script)]
        self.run_command_with_progress(command,cwd=self.dir,prefix='Building JHTagger:',output_width=80)
        # sub.check_call(command,cwd=self.dir,env=env,executable=self.executable,stderr=sub.DEVNULL,stdout=sub.DEVNULL)
        return

    def LoadJHTagger(self,quiet=False):
        # Load our custom ROOT library.
        try:
            a = rt.JHTagger
            return
        except:
            pass

        # Note that we *also* need to fetch some include files -- this has something to do with using the ROOT interpreter.
        # Also note that we have multiple library paths -- to allow for Linux/macOS compatibility.
        build_dir = self.dir + '/jhtagger/build'
        include_dir = self.dir + '/jhtagger/inc'
        custom_lib_paths = [os.path.realpath('{}/lib/libJHTagger.{}'.format(build_dir,x)) for x in self.lib_extensions]
        custom_inc_paths = glob.glob(include_dir + '/**/*.h',recursive=True)

        # Check for any of the libraries.
        found_libary = False
        for libpath in custom_lib_paths:
            if(pathlib.Path(libpath).exists()):
                custom_lib_path = libpath
                found_libary = True
                break

        if(not found_libary):
            if(not quiet): print('Error: The JHTagger lib has not been built!')
            assert False

        for inc_path in custom_inc_paths:
            command = '#include "{}"'.format(inc_path)
            status = rt.gInterpreter.Declare(command)
            assert status, 'The following header file did not load properly: {}'.format(inc_path)

        status = rt.gSystem.Load(custom_lib_path)
        assert status == 0, 'The JHTagger lib did not load properly.'
        return

    def run_command_with_progress(self, command, prefix, cwd=None,
                                show_output_lines=10, output_width=80, output_length=None):
            """
            Run command with progress bar and scrolling output display

            Args:
                command: Command to run
                cwd: Current working directory for command
                prefix: Progress bar prefix
                show_output_lines: Number of recent output lines to show (0 to disable)
                output_width: Width of output display area
                output_length: Optional expected length of output; to use if output has no progress tracking (e.g. "[1%] ...")
            """
            # Pattern to match progress indicators
            progress_pattern = re.compile(r'\[\s*(\d+)%\]')

            # Buffer to store recent output lines
            recent_lines = deque(maxlen=show_output_lines * 2)  # Store more than we show for better filtering

            # Track current progress
            current_progress = 0

            def read_stream(stream, file_handle, line_queue, stream_name):
                """Helper function to read from stdout/stderr in separate threads"""
                try:
                    for line in iter(stream.readline, ''):
                        if line:
                            file_handle.write(line)
                            file_handle.flush()
                            line_queue.put((stream_name, line))
                except Exception as e:
                    line_queue.put(('error', f"Error reading {stream_name}: {str(e)}\n"))
                finally:
                    stream.close()

            with open(self.logfile, 'a') as f, open(self.errfile, 'a') as g:

                # add an empty space -- nice to separate outputs from different commands
                f.write('\n')
                g.write('\n')

                process = sub.Popen(
                    command,
                    cwd=cwd,
                    stdout=sub.PIPE,
                    stderr=sub.PIPE,
                    universal_newlines=True,
                    bufsize=1  # Line buffered
                )

                # Create queues for inter-thread communication
                line_queue = queue.Queue()

                # Start threads to read stdout and stderr
                stdout_thread = threading.Thread(
                    target=read_stream,
                    args=(process.stdout, f, line_queue, 'stdout')
                )
                stderr_thread = threading.Thread(
                    target=read_stream,
                    args=(process.stderr, g, line_queue, 'stderr')
                )

                stdout_thread.daemon = True
                stderr_thread.daemon = True
                stdout_thread.start()
                stderr_thread.start()

                # Initialize display area
                if show_output_lines > 0:
                    # Reserve space for progress bar + output lines
                    for i in range(3 + show_output_lines):
                        print()
                    # Initial display
                    printProgressWithOutput(current_progress, 100, recent_lines,
                                        prefix=prefix, suffix='Complete',
                                        max_lines=show_output_lines, width=output_width)
                else:
                    # Just show progress bar
                    printProgressBar(current_progress, 100, prefix=prefix, suffix='Complete')

                counter = 0
                try:
                    while process.poll() is None or not line_queue.empty():
                        try:
                            # Get line with timeout to check process status
                            stream_name, line = line_queue.get(timeout=0.1)

                            # Add to recent lines buffer
                            recent_lines.append(line)

                            # Update progress -- either checking for progress indicator,
                            # or expecting some fixed length of output.
                            if(output_length is not None):
                                current_progress = (counter / output_length) * 100.

                            # Check for progress indicator in stdout
                            elif stream_name == 'stdout':
                                match = progress_pattern.search(line)
                                if match:
                                    current_progress = int(match.group(1))

                            # Update display (both progress and output)
                            if show_output_lines > 0:
                                printProgressWithOutput(current_progress, 100, recent_lines,
                                                    prefix=prefix, suffix='Complete',
                                                    max_lines=show_output_lines, width=output_width)
                            else:
                                printProgressBar(current_progress, 100, prefix=prefix, suffix='Complete')

                            counter += 1
                        except queue.Empty:
                            continue  # Continue checking if process is still running

                except KeyboardInterrupt:
                    process.terminate()
                    raise

                # Wait for threads to complete
                stdout_thread.join(timeout=1.0)
                stderr_thread.join(timeout=1.0)

                # Final progress update
                if show_output_lines > 0:
                    printProgressWithOutput(100, 100, recent_lines,
                                        prefix=prefix, suffix='Complete',
                                        max_lines=show_output_lines, width=output_width)
                    print()  # Move cursor to next line after final display
                else:
                    printProgressBar(100, 100, prefix=prefix, suffix='Complete')
                    print()  # Move to next line

                # Wait for process to complete and check return code
                return_code = process.wait()
                if return_code != 0:
                    raise sub.CalledProcessError(return_code, command)


class JohnsHopkinsTagger:
    """
    The Johns Hopkins top tagger. This class leverages our custom
    ROOT::JHTagger::JohnnyTagger() class, that interfaces with the Fastjet
    JH tagger implementation.

    See: https://arxiv.org/abs/0806.0848 [Phys.Rev.Lett. 101 (2008) 142001]

    """
    def __init__(self,delta_p=0.1,delta_r=0.19,cos_theta_W_max=0.7,top_mass_range=(150.,200.),W_mass_range=(65.,95.), mode='filter',tag_name=None, n_w_constituents_max=100):
        self.mode = mode
        assert self.mode in ['tag','filter']
        self.tag_name = tag_name
        self.w_name = None
        self.w_nconst_name = None
        self.w_constituents_name = None
        self.w_constituents_name = None
        # Prepare the JHTagger library.
        self.jh_setup = None

        self.tagger = None
        self.SetDeltaP(delta_p,init=False)
        self.SetDeltaR(delta_r,init=False)
        self.SetCosThetaWMax(cos_theta_W_max,init=False)
        self.SetTopMassRange(*top_mass_range,init=False)
        self.SetWMassRange(*W_mass_range,init=False)
        # self.InitTagger()

        # Transient, per-jet variables
        self.tag_status = False
        self.w_candidate = None
        self.w_candidate_cyl = None
        self.w_constituent = None
        self.w_constituent_cyl = None

        # Event-level; one extra dim w.r.t. above
        # TODO: Rework this? Might be confusing naming scheme
        self.tags = None
        self.w_candidates = None
        self.w_candidates_cyl = None
        self.n_constituents = None
        self.w_constituents = None
        self.w_constituents_cyl = None

        self.n_w_constituents_max = n_w_constituents_max

        self.print_prefix = '\n\t\tJohnsHopkinsTagger'
        self.citations = {
            "JohnsHopkinsTopTagger":
            """
@article{Kaplan:2008ie,
    author = "Kaplan, David E. and Rehermann, Keith and Schwartz, Matthew D. and Tweedie, Brock",
    title = "{Top Tagging: A Method for Identifying Boosted Hadronically Decaying Top Quarks}",
    eprint = "0806.0848",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    doi = "10.1103/PhysRevLett.101.142001",
    journal = "Phys. Rev. Lett.",
    volume = "101",
    pages = "142001",
    year = "2008"
}
            """
        }
    def GetCitations(self):
        return self.citations


    def SetDeltaP(self,p,init=True):
        self.delta_p = p
        if(init): self.InitTagger()

    def SetDeltaR(self,r,init=True):
        self.delta_r = r
        if(init): self.InitTagger()

    def SetCosThetaWMax(self,c,init=True):
        self.cos_theta_W_max = c
        if(init): self.InitTagger()

    def SetTopMassRange(self,mmin,mmax,init=True):
        self.top_mass_range = (mmin,mmax)
        if(init): self.InitTagger()

    def SetWMassRange(self,mmin,mmax,init=True):
        self.W_mass_range = (mmin,mmax)
        if(init): self.InitTagger()

    def InitTagger(self):
        if(self.tagger is not None):
            del self.tagger # TODO: Is this necessary?

        self.tagger = rt.JHTagger.JohnnyTagger(self.delta_p,self.delta_r,self.cos_theta_W_max,*self.top_mass_range,*self.W_mass_range)

    def _tag(self,vecs):
        #NOTE: This usage of TagJet() is a bit awkward, but I've had some issues when trying to pass a fastjet.PseudoJet object.
        #      I think this has to do with some weirdness around the Fastjet Python interface, or Python/C++ interfaces in general.
        self.tagger.TagJet(
            rt.std.vector('double')(vecs[:,0].flatten()),
            rt.std.vector('double')(vecs[:,1].flatten()),
            rt.std.vector('double')(vecs[:,2].flatten()),
            rt.std.vector('double')(vecs[:,3].flatten())
        )
        self.tag_status = self.tagger.GetStatus()
        self.w_candidate = None
        if(self.tag_status):
            try:
                w = self.tagger.GetWCandidate() # Fastjet::PseudoJet (C++ type, PyROOT seems to handle interface here!)
                self.w_candidate     = np.array([w.e(),w.px(),w.py(),w.pz()  ]) # TODO: Can this array return be handled by the C++/ROOT class?
                self.w_candidate_cyl = np.array([w.pt(),w.eta(),w.phi(),w.m()])
            except:
                self.w_candidate     = np.full(4,np.nan)
                self.w_candidate_cyl = np.full(4,np.nan)
        return

    def _getWConstituents(self):
        if(not self.tag_status):
            print('BLAM')
            return None
        pt = np.array(self.tagger.GetWCandidateConstituentsProperty("pt"))
        n = pt.shape[0]
        if(n == 0):
            return
        ordering = np.argsort(-pt)
        self.w_constituent     = np.vstack([self.tagger.GetWCandidateConstituentsProperty(x) for x in ["E","px","py","pz"]]  ).T[ordering]
        self.w_constituent_cyl = np.vstack([self.tagger.GetWCandidateConstituentsProperty(x) for x in ["pt","eta","phi","m"]]).T[ordering]

    def ModifyInitialization(self, obj : 'JetFinder'):
        #NOTE: Might want to move this to constructor, which will need to take obj as input.
        self.jh_setup = JHTaggerSetup(obj.configurator)
        # self.jh_setup.SetConfigurator(obj.configurator)
        self.jh_setup.FullPreparation()
        self.InitTagger()
        return

    def ModifyInputs(self,obj : 'JetFinder'):
        return

    def ModifyJets(self, obj : 'JetFinder'):
        """
        This function will tag jets with the JH tagger, and fill the corresponding branches.
        """
        import fastjet as fj # NOTE: In practice, fastjet will have been initialized already by JetFinder. Can similarly do this in Softdrop

        # Fetch the jet constituents, just to be safe -- this makes sure that they are up-to-date.
        obj._fetchJetConstituents()

        self.tags = {i:False for i in obj.jets_dict.keys()}
        self.w_candidates = {i:np.zeros(4) for i in obj.jets_dict.keys()}
        self.w_candidates_cyl = {i:np.zeros(4) for i in obj.jets_dict.keys()}
        self.n_constituents = {i:0 for i in obj.jets_dict.keys()}
        self.w_constituents = {i:np.zeros((self.n_w_constituents_max,4)) for i in obj.jets_dict.keys()}
        self.w_constituents_cyl = {i:np.zeros((self.n_w_constituents_max,4)) for i in obj.jets_dict.keys()}

        for key in obj.jets_dict.keys():
            self._tag(obj.constituent_vectors[key]) # fills self.tag_status, self.w_candidate
            self.tags[key] = self.tag_status
            if(self.tag_status): # if not tagged, there is no W -- so we can safely skip filling
                self._getWConstituents() # fills self.w_constituent, self.w_constituent_cyl
                self.n_constituents[key] = len(self.w_constituent)
                embed_array_inplace(self.w_constituent,self.w_constituents[key])
                embed_array_inplace(self.w_constituent_cyl,self.w_constituents_cyl[key])

        if(self.mode=='filter'):
            obj.jet_ordering = [key for key in obj.jet_ordering if self.tags[key]]
            obj._updateJetDictionary()
            # Refresh vectors and constituents -- always need to do this if we filter jets_dict.
            obj._jetsToVectors()
            obj._fetchJetConstituents()

    def ModifyConstituents(self, obj : 'JetFinder'):
        return

    def ModifyWrite(self,obj : 'JetFinder'):
        if(self.mode=='filter'):
            return # do nothing
        else:
            self._initializeBuffer(obj) # will initialize buffer if it doesn't already exist
            self._addFlagToBuffer(obj)
            self._addWToBuffer(obj)

    def _initializeBuffer(self,obj : 'JetFinder'):
        """
        Used if self.mode=='tag', in which case we're writing all jets,
        and including a new branch to indicate whether or not a jet is
        JH-tagged.
        """
        self._createBranchNames(obj)

        if(self.tag_name not in obj.buffer.keys()):
            obj.buffer.create_array(self.tag_name,(obj.n_jets_max,),dtype=bool)

        if(self.w_name not in obj.buffer.keys()):
            obj.buffer.create_array(self.w_name,         (obj.n_jets_max,4),dtype=np.dtype('f8'))
            obj.buffer.create_array(self.w_name + '_cyl',(obj.n_jets_max,4),dtype=np.dtype('f8'))

        if(self.w_constituents_name not in obj.buffer.keys()):
            obj.buffer.create_array(self.w_nconst_name, (obj.n_jets_max),dtype=np.dtype('i4'))
            obj.buffer.create_array(self.w_constituents_name,          (obj.n_jets_max,self.n_w_constituents_max,4),dtype=np.dtype('f8'))
            obj.buffer.create_array(self.w_constituents_name + '_cyl', (obj.n_jets_max,self.n_w_constituents_max,4),dtype=np.dtype('f8'))
        return

    def _createBranchNames(self,obj : 'JetFinder'):
        if(self.tag_name is None):
            self.tag_name = '{}.JHTag'.format(obj.jet_name)

        # Also create keys corresponding to the W-boson candidate, and its constituents,
        # that are identified by the JH tagger.
        self.w_name = '{}.WBosonCandidate.Pmu'.format(self.tag_name)
        self.w_nconst_name = '{}.WBosonCandidate.Constituents.N'.format(self.tag_name)
        self.w_constituents_name = '{}.WBosonCandidate.Constituents.Pmu'.format(self.tag_name)

    def _addFlagToBuffer(self,obj : 'JetFinder'):
        """
        Adds the JH tags to the buffer, for writing.
        Note that the pT sorting of obj is applied,
        which will have been filled by obj._ptSort().
        """
        obj.buffer.set(self.tag_name,obj._i,[self.tags[i] for i in obj.jet_ordering])

    def _addWToBuffer(self,obj : 'JetFinder'):
        """
        Adds the JH W candidate info to the buffer.
        Note that the pT sorting of obj is applied,
        which will have been filled by obj._ptSort().
        """
        obj.buffer.set(self.w_name,obj._i,[self.w_candidates[i] for i in obj.jet_ordering])
        obj.buffer.set(self.w_name + '_cyl',obj._i,[self.w_candidates_cyl[i] for i in obj.jet_ordering])

        obj.buffer.set(self.w_nconst_name,obj._i,[self.n_constituents[i] for i in obj.jet_ordering])
        obj.buffer.set(self.w_constituents_name,obj._i,[self.w_constituents[i] for i in obj.jet_ordering])
        obj.buffer.set(self.w_constituents_name + '_cyl',obj._i,[self.w_constituents_cyl[i] for i in obj.jet_ordering])

    def _print(self,val):
        print('{}: {}'.format(self.print_prefix,val))
        return
