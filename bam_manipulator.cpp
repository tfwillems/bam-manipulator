
#include <getopt.h>

#include <cstdlib>
#include <iostream>
#include <set>
#include <sstream>
#include <queue>
#include <unistd.h>

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
	    << "\t" << "--bam     <input_bam_file> "  << "\t" << "Input BAM file path"               << "\n"
	    << "\t" << "--out     <output_bam_file>"  << "\t" << "Output BAM file path"              << "\n"

	    << "Optional parameters:" << "\n"
	    << "\t" << "--help                     "  << "\t" << "Print this help message and exit"  << "\n"
	    << "\t" << "--version                  "  << "\t" << "Print HipSTR version and exit"     << "\n"
	    << std::endl;
}

void parse_command_line_args(int argc, char** argv, std::string& input_bam_file, std::string& output_bam_file){
  if (argc == 1 || (argc == 2 && std::string("-h").compare(std::string(argv[1])) == 0)){
    print_usage();
    exit(0);
  }

  int print_help = 0, print_version = 0;
  static struct option long_options[] = {
    {"h",       no_argument, &print_help,    1},
    {"help",    no_argument, &print_help,    1},
    {"version", no_argument, &print_version, 1},
    {"bam",     required_argument, 0, 'b'},
    {"out",     required_argument, 0, 'o'},
    {0, 0, 0, 0}
  };

  int c;
  while (true){
    int option_index = 0;
    c = getopt_long(argc, argv, "b:o:", long_options, &option_index);
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
      input_bam_file  = std::string(optarg);
      break;
    case 'o':
      output_bam_file = std::string(optarg);
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
  std::string input_bam_file = "", output_bam_file = "";
  parse_command_line_args(argc, argv, input_bam_file, output_bam_file);

  if (input_bam_file.empty())
    printErrorAndDie("--bam option required");
  if (output_bam_file.empty())
    printErrorAndDie("--out option required");

  BamTools::BamReader bam_reader;
  if (!bam_reader.Open(input_bam_file))
    printErrorAndDie("Failed to open the input BAM file");

  BamTools::BamWriter bam_writer;
  BamTools::RefVector ref_vector = bam_reader.GetReferenceData();
  if (!bam_writer.Open(output_bam_file, bam_reader.GetHeaderText(), ref_vector))
    printErrorAndDie("Failed to open the output BAM file");

  bam_reader.Rewind();
  std::set<std::pair<std::string, bool> > processed_ids;
  BamTools::BamAlignment alignment;
  int32_t filt_count = 0, pass_count = 0;
  while (bam_reader.GetNextAlignment(alignment)){
    std::pair<std::string, bool> key(alignment.Name, alignment.IsFirstMate());
    if (processed_ids.find(key) != processed_ids.end()){
      filt_count++;
      continue;
    }
    else {
      pass_count++;
      bam_writer.SaveAlignment(alignment);
      processed_ids.insert(key);
    }
  }

  std::cerr << "Filtered out " << filt_count << " out of " << (filt_count+pass_count) << " reads due to multiple split alignments" << std::endl;
  bam_reader.Close();
  bam_writer.Close();
  return 0;
}
