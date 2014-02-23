#include "api/BamMultiReader.h"
#include "api/BamWriter.h"
#include <iostream>
using namespace BamTools;
using namespace std;

int main() 
{
  // at some point, start our merge operation
  vector<string> inputFilenames;

  inputFilenames.push_back("test.bam");
  string outputFilename = "test.out";
  // provide some input & output filenames
  // attempt to open our BamMultiReader
  BamMultiReader reader;
  if ( !reader.Open(inputFilenames) ) {
      cerr << "Could not open input BAM files." << endl;
      return 1;
  }
  // retrieve 'metadata' from BAM files, these are required by BamWriter
  const SamHeader header = reader.GetHeader();
  const RefVector references = reader.GetReferenceData();
  // attempt to open our BamWriter
  BamWriter writer;
  if ( !writer.Open(outputFilename, header, references) ) {
      cerr << "Could not open output BAM file" << endl;
      return 1;
  }
  // iterate through all alignments, only keeping ones with high map quality
  BamAlignment al;
  while ( reader.GetNextAlignmentCore(al) ) {
      if ( al.MapQuality >= 90 )
	  writer.SaveAlignment(al);
      }
  // close the reader & writer
  reader.Close();
  writer.Close();
  
  for (int i = 0; i < inputFilenames.size(); i++) {
      cout << "Opening " << endl;
      cerr << inputFilenames[i] << endl;
  }
  return 0;
}

// // merge is now complete, continue whatever we were doing
