#include "Graph.h"
#include "Utility.h"
#include "Timer.h"

using namespace std;
#define LEN_LIMIT (1<<10)
char filename[LEN_LIMIT];
int k, maxsec = 1800;

signed main(int argc, char *argv[]) {
	printf("\n-----------------------------------------------------------------------------------------\n");
	if (argc == 3) {
		strncpy(filename, argv[1], LEN_LIMIT);
		k = atoi(argv[2]);
		Graph *graph = new Graph(filename, k);
		graph->read();
		graph->search();
		graph->write();
		// delete graph; // there are some bugs in releasing memory
	}
	else printf("[usage]: binary filename k\n");
	printf("-----------------------------------------------------------------------------------------\n\n");
}