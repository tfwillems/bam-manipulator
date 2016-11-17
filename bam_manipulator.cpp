#include <getopt.h>

#include <cstdlib>
#include <iostream>
#include <set>
#include <sstream>
#include <queue>

#include "error.h"
#include "stringops.h"
#include "version.h"

#include "bamtools/include/api/BamAlignment.h"
#include "bamtools/include/api/BamMultiReader.h"
#include "bamtools/include/api/BamWriter.h"


bool file_exists(std::string path){
  return (access(path.c_str(), F_OK) != -1);
}

void print_usage(){
  std::cerr << "Usage: BamManipulator --bam <bam_file>" << "\n" << "\n"

	    << "Required parameters:" << "\n"
	    << "\t" << "--bam     <bam_file>       "  << "\t" << "BAM file path"                              << "\n"

	    << "Optional parameters:" << "\n"
	    << "\t" << "--help                    "  << "\t" << "Print this help message and exit"            << "\n"
	    << "\t" << "--version                 "  << "\t" << "Print HipSTR version and exit"               << "\n"
	    << std::endl;
}

void parse_command_line_args(int argc, char** argv, std::string& bam_file_string){
  if (argc == 1 || (argc == 2 && std::string("-h").compare(std::string(argv[1])) == 0)){
    print_usage();
    exit(0);
  }

  int print_help = 0, print_version = 0;
  static struct option long_options[] = {
    {"h",       no_argument, &print_help,    1},
    {"help",    no_argument, &print_help,    1},
    {"version", no_argument, &print_version, 1},
    {"bams",    required_argument, 0, 'b'},
    {0, 0, 0, 0}
  };

  int c;
  while (true){
    int option_index = 0;
    c = getopt_long(argc, argv, "b:", long_options, &option_index);
    if (c == -1)
      break;

    if (optarg != NULL){
      std::string val(optarg);
      if (string_starts_with(val, "--"))
	printErrorAndDie("Argument to option --" + std::string(long_options[option_index].name) + " cannot begin with \"--\"\n\tBad argument: " + val);
    }

    switch(c){
    case 0:
      break;
    case 'b':
      bam_file_string = std::string(optarg);
      break;
    case '?':
      exit(1);
      break;
    default:
      abort();
      break;
    }
  }

  if (optind < argc) {
    std::stringstream msg;
    msg << "Did not recognize the following command line arguments:" << "\n";
    while (optind < argc)
      msg << "\t" << argv[optind++] << "\n";
    msg << "Please check your command line syntax or type ./BamManipulator --help for additional information" << "\n";
    printErrorAndDie(msg.str());
  }

  if (print_version == 1){
    std::cerr << "BamManipulator version " << VERSION << std::endl;
    exit(0);
  }

  if (print_help){
    print_usage();
    exit(0);
  }
}



int main(int argc, char** argv){
  std::string bam_file_string = "";
  parse_command_line_args(argc, argv, bam_file_string);

  if (bam_file_string.empty())
    printErrorAndDie("--bam option required");

  BamTools::BamReader bam_reader;
  if (!bam_reader.Open(bam_file_string))
    printErrorAndDie("Failed to the input BAM file");

  // TO DO: Check that the BAM is sorted by read ID


  bam_reader.Rewind();
  std::vector<BamTools::BamAlignment> alignments;
  BamTools::BamAlignment alignment;
  int32_t filt_count = 0, pass_count = 0;
  std::string prev_id = "";
  while(bam_reader.GetNextAlignment(alignment)){
    if (alignment.Name.compare(prev_id) != 0){
      if (!alignments.empty()){
	if (alignments.size() == 2){
	  pass_count += alignments.size();
	}
	else
	  filt_count += alignments.size();
      }

      alignments.clear();
      alignments.push_back(alignment);
      prev_id = alignment.Name;
    }
    else
      alignments.push_back(alignment);
  }

  if (alignments.size() == 2){
    pass_count += alignments.size();
  }
  else
    filt_count += alignments.size();

  std::cerr << "Filtered out " << filt_count << " out of " << (filt_count+pass_count) << " alignments" << std::endl;
  bam_reader.Close();
  return 0;
}
