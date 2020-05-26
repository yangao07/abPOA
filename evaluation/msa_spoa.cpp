/* To complie:
   g++ -O3 msa_spoa.cpp -std=c++11 -I include/ -L build/lib/ -lspoa -o msa_spoa
  */

#include <getopt.h>
#include <stdio.h>
#include <cstdint>
#include <iostream>
#include <fstream>
#include <string>
#include <sys/time.h>
#include <time.h>
#include <assert.h>
#include <exception>
#include <getopt.h>
#include "spoa/spoa.hpp"
#include "spoa/graph.hpp"
#include "spoa/alignment_engine.hpp"

using namespace std;
using Alignment = std::vector<std::pair<std::int32_t, std::int32_t>>;

// read fasta format
string get_seq(ifstream& in_file)
{
    char c = in_file.get();
    if (in_file.eof()) return "";
    else if (c != '>') {
        cout << "Unexpected format." << endl;
        exit(1);
    }
    string seq;
    getline(in_file, seq); // header line
    getline(in_file, seq); // sequence line
    return seq;
}

string get_seq_for_racon(ifstream& in_file, bool *start)
{
    char c = in_file.get();
    if (in_file.eof()) return "";
    else if (c != '>') {
        cout << "Unexpected format." << endl;
        exit(1);
    }
    string seq;
    getline(in_file, seq); // header line
    if (seq[0] == '0') *start = true;
    else *start = false;
    getline(in_file, seq); // sequence line
    return seq;
}

void help() {
    std::cout <<
        "\n"
        "usage: ./msa_spoa [options ...]\n"
        "\n"
        "    options:\n"
        "        -m <int>\n"
        "            default: 2\n"
        "            score for matching bases\n"
        "        -x <int>\n"
        "            default: 4\n"
        "            penalty for mismatching bases\n"
        "        -o <int(,int)>\n"
        "            default: gap_open1 = 4, gap_open2 = 24\n"
        "            gap opening penalty (must be non-negative)\n"
        "        -e <int(,int)>\n"
        "            default: gap_ext1 = 2, gap_ext2 = 1\n"
        "            gap extension penalty (must be non-negative)\n"
        "        -s <file>\n"
        "            default: seq.fa\n"
        "            the input sequence set\n"
        "        -n <int>\n"
        "            default: 10\n"
        "            number of sequences in each set\n"
        "        -r \n"
        "            input sequence set is from racon\n"
        "        -h \n"
        "            prints the usage\n";
}

int main(int argc, char** argv) {
    string seq_file = "seq.fa";

    std::uint8_t algorithm = 1;
    std::int8_t m = 2;
    std::int8_t x = -4;
    std::int8_t o1 = -4;
    std::int8_t e1 = -2;
    std::int8_t o2 = -24;
    std::int8_t e2 = -1;

    char opt, *s; int n_seqs = 0, for_racon = 0;
    while ((opt = getopt(argc, argv, "l:m:x:o:n:e:q:c:s:rh")) != -1) {
        switch (opt) {
            case 'm': m = atoi(optarg); break;
            case 'x': x = 0-atoi(optarg); break;
            case 'o': o1 = 0-strtol(optarg, &s, 10); if (*s == ',') o2 = 0-strtol(s+1, &s, 10); break;
            case 'e': e1 = 0-strtol(optarg, &s, 10); if (*s == ',') e2 = 0-strtol(s+1, &s, 10); break;
            case 'n': n_seqs = atoi(optarg); break;
            case 's': seq_file = optarg; break;
            case 'r': for_racon = 1; break;
            case 'h': help(); return 0;
            default: help(); return 1;
        }
    }

    std::int8_t oe1=o1+e1, oe2=o2+e2;

    std::unique_ptr<spoa::AlignmentEngine> alignment_engine;
    try {
        alignment_engine = spoa::createAlignmentEngine(
                static_cast<spoa::AlignmentType>(algorithm), m, x, oe1, e1, oe2, e2);
    } catch(std::invalid_argument& exception) {
        std::cerr << exception.what() << std::endl;
        return 1;
    }

    ifstream fp_seq;
    fp_seq.open(seq_file, ios::in);
    assert(fp_seq.is_open());

    struct timeval start_time, end_time;
    double runtime = 0; int seq_i; string seq;

    if (for_racon) {
        seq_i = 0;
        int init = 1; bool start;
        while(1) {
            auto graph = spoa::createGraph();

            while(1) {
                seq = get_seq_for_racon(fp_seq, &start);
                if(seq == "" || start) break;

                gettimeofday(&start_time, NULL);
                auto alignment = alignment_engine->align(seq, graph);
                graph->add_alignment(alignment, seq);
                gettimeofday(&end_time, NULL);
                runtime += (end_time.tv_sec - start_time.tv_sec) + (end_time.tv_usec - start_time.tv_usec) * 1e-6;
                seq_i++;
            }
            if (init) init = 0;
            else {
                // cout << seq_i << endl;

                gettimeofday(&start_time, NULL);
                std::string consensus = graph->generate_consensus();
                std::cout << ">Consensus_sequence" << endl;
                std::cout << consensus.c_str() << endl;
                gettimeofday(&end_time, NULL);
                runtime += (end_time.tv_sec - start_time.tv_sec) + (end_time.tv_usec - start_time.tv_usec) * 1e-6;
                if (seq_i == 0 || fp_seq.eof()) break;
            }
            if (start) {
                gettimeofday(&start_time, NULL);
                auto alignment = alignment_engine->align(seq, graph);
                graph->add_alignment(alignment, seq);
                gettimeofday(&end_time, NULL);
                runtime += (end_time.tv_sec - start_time.tv_sec) + (end_time.tv_usec - start_time.tv_usec) * 1e-6;
                seq_i = 1;
            }
        }
    } else {
        while(1) {
            auto graph = spoa::createGraph();

            seq_i = 0;
            while(1) {
                seq = get_seq(fp_seq);
                if(seq == "") break;

                gettimeofday(&start_time, NULL);
                auto alignment = alignment_engine->align(seq, graph);
                graph->add_alignment(alignment, seq);
                gettimeofday(&end_time, NULL);
                runtime += (end_time.tv_sec - start_time.tv_sec) + (end_time.tv_usec - start_time.tv_usec) * 1e-6;

                if (++seq_i == n_seqs) break;
            }

            if (seq_i == 0 || fp_seq.eof()) break;

            gettimeofday(&start_time, NULL);
            std::string consensus = graph->generate_consensus();
            std::cout << ">Consensus_sequence" << endl;
            std::cout << consensus.c_str() << endl;
            gettimeofday(&end_time, NULL);
            runtime += (end_time.tv_sec - start_time.tv_sec) + (end_time.tv_usec - start_time.tv_usec) * 1e-6;
        }
    }
    fprintf(stderr, "%.2f ", runtime);

    fp_seq.close();
    return 0;
}
