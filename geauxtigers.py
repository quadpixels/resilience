# This is module geauxtigers.
# It contains some utilities for injection.

import os, subprocess, threading, signal

def Source(script, update=1):
	pipe = subprocess.Popen(". %s; env" % script, stdout=subprocess.PIPE, shell=True)
	data = pipe.communicate()[0];
	lines=[]
	for line in data.splitlines():
		l=line.split("=", 1)
		if(len(l)==2): 
			lines.append(l)
	env = dict(lines)
	print "%d environment variables." % len(lines)
	if update:
		os.environ.update(env)
	return env

def runTimed(timeout, cmd):
	event = threading.Event()
	pipe = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=False)
	killer = KillerThread(pipe.pid, pipe, timeout, event)
	killer.start()
	(out, err) = pipe.communicate()
	event.set()
	return (out, err, pipe.returncode)

class KillerThread(threading.Thread):
	def __init__(self, pid, pipe, timeout, event):
		threading.Thread.__init__(self)
		self.pid = pid
		self.timeout = timeout
		self.event = event
		self.setDaemon(True)
		self.pipe = pipe
	
	def run(self):
		self.event.wait(self.timeout)
		if not self.event.isSet():
			try:
				os.kill(self.pid, signal.SIGKILL)
				self.pipe.kill();
				print "Killed -- ran out of time!"
			except OSError, e:
				print "OSError"
				pass

def writeHistogram(array, filename):
	fh = open(filename, 'w')
	for a in array:
		fh.write(str(a) + "\n")
	fh.close()

# Definition of some regexes
re_poecc_encd = r">> Elems encoded: (?P<elems_encd>[0-9]+).*";
re_poecc_corrd = r">> Elems corrected: (?P<elems_corrd>[0-9]+).*";
re_myint_1 = "\[MyInt](.*\t)+.*";
re_myptr_1 = "\[MyPtr](.*\t)+.*";
re_dump = ">>.*";
re_entered_handler = r" >> Caught SIGSEGV signal \((?P<sigsegv_num>[0-9]+) out of [0-9]+ allowed\)"
re_poke_instcnt = "__POKE__ InstCount: (?P<instcnt>[0-9]+).*";
re_poke_rndnums = "__POKE__ Generated (?P<rndnums>[0-9]+) random numbers";
re_sig_summary = "\[MySigsegvHandler\].*"                        
re_fir = r">> FastLowPassFir correct \((?P<fir_correct>[0-9]+)/(?P<fir_incorrect>[0-9]+)/(?P<fir_total>[0-9]+)/(?P<fir_recalc>[0-9]+)\).*";
re_fftrecalc = r">> FFT recalculation: (?P<fft_recalc>[0-9]+).*";
re_pass4 = r">> Pass 4 recalculation: (?P<pass4_recalc>[0-9]+).*";
re_stopwatch = r"\[PID=(?P<pid>[0-9]+)\]\[my_stopwatch_show #(?P<sw_id>[0-9]+)\] (?P<microsec>[0-9]+) microseconds \((?P<desc>.*)\)\.";
re_num_retries = "\[FT3\] # of retries - (?P<routine>.*): (?P<retries>[0-9]+).*"
re_ft3test_rmsd = "\[FT3Test (?P<routine>.*)\].*RMSD=(?P<rmsd>.*)"
