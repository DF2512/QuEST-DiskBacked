# QuEST-DiskBacked

Preliminary setup may be required such as installing liburing (sudo apt install liburing-dev)

## Testing guide
- **Edit the files "example/tutorials/min_example.cpp" (line 48) and "examples/tutorials/min_example1.cpp" (line 109) to add target directory path for storing state vector amplitudes as a string**
- If you encounter the error 'Floating point exception (core dumped)', please refer to the above point
- Follow QuEST guidance for configuring and building
- Run min_example to execute a simple runtime benchmark
- Run min_example1 to verify accuracy of implementation

## Implementation Overview

No QuEST source files were altered other than quest.h,in (adding new headers to the main header), and CMakeLists.txt in order to add liburing as a library and ensure min_example executables were generated.

### New original files
- quest/src/core/diskbackedstate.cpp (include/.h): contains metadata of state vector, functions for I/O, and functions for initialisation/probability/measurement
- quest/src/core/gatescheduler.cpp (include/.h): contains the gate schedule object, functions that call the QuEST gate functions, and the permutation/partition generation
- quest/src/core/chunkmanager.cpp (include/.h): contains transition information generation and stores metadata regarding the current permutation and chunkmap
- quest/src/core/run.cpp (include/.h): contains the triple buffered pipeline execution function
- quest/src/core/hook.cpp (include/.h): contains the temporary Qureg function, enabling the implementation to interact with QuEST functions
- min_example.cpp and min_example1.cpp: as above (note here than min_example.c is NOT original)

## Usage Guide
First, initialise the QuEST environment as usual. Also include required libraries and files. (TODO: figure out which, if any, of these includes are not needed. chrono is only needed when timing for benchmarking purposes as an example).
```cpp
#include "quest.h"
#include "diskbackedstate.h"
#include "gatescheduler.h"
#include "run.h"
#include "chunkmanager.h"
#include <vector>
#include <string>
#include <iostream>
#include <cmath>
#include <random>
#include <ctime>
#include <complex>
#include <algorithm>
#include <iomanip>
#include <numeric>
#include <exception>
#include <omp.h>
#include <set>
#include <limits>
#include <chrono>
#include <fstream>
#include <filesystem>

initQuESTEnv();
```
Create the object that handles the disk and the gate schedule. Also create the vector that holds the directory path(s). The parameters for the object are as follows:
- numQubits: the total number of qubits.
- numBlocks: the number of blocks to split the state vector into.
- chunksPerBlock: the number of chunks each block is split into.
- roots: the vector of directory strings.
- maxBlocksInMemory: the maximum number of blocks allowed in memory at any one time.
```cpp
std::vector<std::string> roots = {
	"directory1/.../.../path",
	"directory2/.../.../path"
}

DiskBackedState state(numQubits, numBlocks, chunksPerBlock, roots, maxBlocksInMemory);

GateScheduler schedule;
```
Be careful when selecting the values for each parameter. There are some important guidelines to follow in order to ensure that no unexpected errors occur. 

First, numBlocks and chunksPerBlock should both be powers of 2 (i.e. 1, 2, 4, 8,...) to ensure that each chunk and each block contains a whole number of qubits. 

Next, it is important to calculate the storage size of each block and chunk in order to adhere to restrictions imposed by liburing. It is simple to calculate these sizes, and can be represented both in bytes and in the number of qubits. In bytes, divide the size of the full state vector in bytes by the number of blocks to get the size of a block in bytes, and divide the size of a block by the number of chunks per block to get the size of a chunk. To represent in qubits, subtract the log2(numBlocks) from the numQubits to get the number of qubits represented by each block, and subtract the log2(chunksPerBlock) from the number of qubits per block to get the number of qubits represented by each chunk.

The limitations imposed are as follows:
- Ensure that the size in bytes of a block multiplied by the maxBlocksInMemory is less than or equal to the available RAM.
- Ensure that the size in bytes of each chunk is no less than 4,096B, equivalent to 12 qubits (this does mean that state vectors with less than 12 qubits cannot be simulated using this method currently). liburing's registered buffers with O_DIRECT require buffer alignment of 4096, which is guaranteed using the above restriction.
- Ensure that the size in bytes of each chunk is no more than 1GiB, equivalent to 26 qubits. liburing's fixed files must be a maximum of 1GiB, and each chunk is mapped to one file.
- It is recommended that numBlocks and chunksPerBlock are roughly equal to ensure the minimum number of subcircuits are partitioned, maximising performance in theory. If being equal is difficult, it is better to have more chunks per block than number of blocks.
- It is recommended that there are a power of 2 directories, as there will always be a power of 2 number of files, and files are evenly split across the supplied directories. 


Now, the state vector can be initalised, and subsequently passed into a circuit by scheduling the circuit and running the schedule on the state vector.
```cpp
state.diskBacked_initRandomPureState(); // init function names are analogous to base QuEST. Supported
										// init functions are random, zero, plus, and blank.

schedule.addHadmard(target); // Gate function names and params are also analagous to base QuEST.
schedule.addFullQFT(); // The full list of supported gate functions can be found in GateScheduler.cpp/.h.

runCircuit(schedule, state, false); // the false param is legacy code, disabling verbose. TODO remove.
```

After this, diagnostics can be performed such as measurement of a qubit or probability calculation.
```cpp
int outcome = state.diskBacked_applyQubitMeasurement(target);
std::cout << "Outcome of qubit " << target << "measurement: " << outcome << std::endl;

qreal prob = state.diskBacked_calcTotalProbability();
std::cout << "Total probability of state vector: " << prob << std::endl;
```

Please note that initialisation functions and the probability function still require optimisation and are relatively slow. The measurement function is _extremely_ slow but it is very hard to implement an optimised measurement function in this implementation.

The full example script is as follows:
```cpp
#include "quest.h"
#include "diskbackedstate.h"
#include "gatescheduler.h"
#include "run.h"
#include "chunkmanager.h"
#include <vector>
#include <string>
#include <iostream>
#include <cmath>
#include <random>
#include <ctime>
#include <complex>
#include <algorithm>
#include <iomanip>
#include <numeric>
#include <exception>
#include <omp.h>
#include <set>
#include <limits>
#include <chrono>
#include <fstream>
#include <filesystem>

int main() {
	initQuESTEnv();

	std::vector<std::string> roots = {
		"directory1/.../.../path",
		"directory2/.../.../path"
	}

	DiskBackedState state(numQubits, numBlocks, chunksPerBlock, roots, maxBlocksInMemory);

	GateScheduler schedule;

	state.diskBacked_initRandomPureState(); // init function names are analogous to base QuEST. Supported
											// init functions are random, zero, plus, and blank.

	schedule.addHadmard(target); // Gate function names and params are also analagous to base QuEST.
	schedule.addFullQFT(); // The full list of supported gate functions can be found in GateScheduler.cpp/.h.

	runCircuit(schedule, state, false); // the false param is legacy code, disabling verbose. TODO remove.

	int outcome = state.diskBacked_applyQubitMeasurement(target);
	std::cout << "Outcome of qubit " << target << "measurement: " << outcome << std::endl;

	qreal prob = state.diskBacked_calcTotalProbability();
	std::cout << "Total probability of state vector: " << prob << std::endl;

	finalizeQuESTEnv();
	return 0;
}