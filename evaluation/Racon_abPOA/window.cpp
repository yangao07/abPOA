/*!
 * @file window.cpp
 *
 * @brief Window class source file
 */

#include <algorithm>
#include <iostream>

#include "window.hpp"

// AaCcGgTtNn ==> 0,1,2,3,4
unsigned char nst_nt4_table[256] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 5 /*'-'*/, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

namespace racon {

std::shared_ptr<Window> createWindow(uint64_t id, uint32_t rank, WindowType type,
    const char* backbone, uint32_t backbone_length, const char* quality,
    uint32_t quality_length) {

    if (backbone_length == 0 || backbone_length != quality_length) {
        fprintf(stderr, "[racon::createWindow] error: "
            "empty backbone sequence/unequal quality length!\n");
        exit(1);
    }

    return std::shared_ptr<Window>(new Window(id, rank, type, backbone,
        backbone_length, quality, quality_length));
}

Window::Window(uint64_t id, uint32_t rank, WindowType type, const char* backbone,
    uint32_t backbone_length, const char* quality, uint32_t quality_length)
        : id_(id), rank_(rank), type_(type), consensus_(), sequences_(),
        qualities_(), positions_() {

    sequences_.emplace_back(backbone, backbone_length);
    qualities_.emplace_back(quality, quality_length);
    positions_.emplace_back(0, 0);
}

Window::~Window() {
}

void Window::add_layer(const char* sequence, uint32_t sequence_length,
    const char* quality, uint32_t quality_length, uint32_t begin, uint32_t end) {

    if (sequence_length == 0 || begin == end) {
        return;
    }

    if (quality != nullptr && sequence_length != quality_length) {
        fprintf(stderr, "[racon::Window::add_layer] error: "
            "unequal quality size!\n");
        exit(1);
    }
    if (begin >= end || begin > sequences_.front().second || end > sequences_.front().second) {
        fprintf(stderr, "[racon::Window::add_layer] error: "
            "layer begin and end positions are invalid!\n");
        exit(1);
    }

    sequences_.emplace_back(sequence, sequence_length);
    qualities_.emplace_back(quality, quality_length);
    positions_.emplace_back(begin, end);
}

bool Window::generate_consensus(abpoa_para_t* abpt,
    bool trim) {

    if (sequences_.size() < 3) {
        consensus_ = std::string(sequences_.front().first, sequences_.front().second);
        return false;
    }
    std::vector<uint32_t> rank;
    rank.reserve(sequences_.size());
    for (uint32_t i = 0; i < sequences_.size(); ++i) {
        rank.emplace_back(i);
    }

    std::sort(rank.begin() + 1, rank.end(), [&](uint32_t lhs, uint32_t rhs) {
        return positions_[lhs].first < positions_[rhs].first; });
    
    int32_t n_seqs = sequences_.size();
    int *seq_lens = (int*)malloc(sizeof(int) * sequences_.size());
    uint8_t **bseqs = (uint8_t**)malloc(sizeof(uint8_t*) * sequences_.size());
    int32_t i, j; int c;
    for (i = 0; i < n_seqs; ++i) {
        seq_lens[i] = sequences_[i].second;
        bseqs[i] = (uint8_t*)malloc(sizeof(uint8_t) * sequences_[i].second);
        // fprintf(stderr, ">%d, %d\n", i, sequences_[i].second);
        for (j = 0; j < sequences_[i].second; ++j) {
            bseqs[i][j] = nst_nt4_table[sequences_[i].first[j]];
            // fprintf(stderr, "%c", "ACGT"[bseqs[i][j]]);
        }
        // fprintf(stderr, "\n");
    }
    // perform abPOA
    uint8_t **cons_seq; int32_t **cons_cov; int *cons_l, cons_n=0;
    abpoa_t *ab = abpoa_init(); 
    abpoa_res_t res;
    uint32_t offset = 0.01 * sequences_.front().second;
    for (j = 0; j < n_seqs; ++j) { // positions_.first/second: 0-based inc_start/end
        i = rank[j];
        res.graph_cigar = 0, res.n_cigar = 0;
        // fprintf(stderr, "%d: %d %d\n", i, positions_[i].first, positions_[i].second); 
        if (positions_[i].first < offset && positions_[i].second > sequences_.front().second - offset) {
            abpoa_align_sequence_to_graph(ab, abpt, bseqs[i], seq_lens[i], &res);
            abpoa_add_graph_alignment(ab, abpt, bseqs[i], seq_lens[i], res, i, n_seqs);
        } else {
            int32_t inc_beg = positions_[i].first + 2, inc_end = positions_[i].second + 2;
            abpoa_align_sequence_to_subgraph(ab, abpt, inc_beg, inc_end, bseqs[i], seq_lens[i], &res);
            int32_t exc_beg, exc_end;
            if (i != 0) abpoa_subgraph_nodes(ab, inc_beg, inc_end, &exc_beg, &exc_end);
            else exc_beg = 0, exc_end = 0;
            abpoa_add_subgraph_alignment(ab, abpt, exc_beg, exc_end, bseqs[i], seq_lens[i], res, i, n_seqs);
        }
        if (res.n_cigar) free(res.graph_cigar);
    }
    abpoa_generate_consensus(ab, abpt, n_seqs, NULL, &cons_seq, &cons_cov, &cons_l, &cons_n);

    if (cons_n > 0) {
        // extract consensus sequence
        for (j = 0; j < cons_l[0]; ++j) consensus_ = consensus_ + "ACGTN"[cons_seq[0][j]];

        if (type_ == WindowType::kTGS && trim) {
            int begin = 0, end = cons_l[0]-1;
            int32_t average_coverage = (n_seqs - 1) / 2;
            for (; begin < cons_l[0]; ++begin) {
                // printf("BEG: %d %d\n", cons_cov[0][begin], average_coverage);
                if (cons_cov[0][begin] >= average_coverage) break;
            }
            for (; end >= 0; --end) {
                // printf("END: %d %d\n", cons_cov[0][end], average_coverage);
                if (cons_cov[0][end] >= average_coverage) break;
            }
            if (begin >= end) {
                fprintf(stderr, "[racon::Window::generate_consensus] warning: "
                        "contig %lu might be chimeric in window %u!\n", id_, rank_);
            } else {
                consensus_ = consensus_.substr(begin, end - begin + 1);
            }
        }

        for (i = 0; i < cons_n; ++i) {
            free(cons_seq[i]); free(cons_cov[i]);
        } free(cons_seq); free(cons_cov); free(cons_l);
        abpoa_free(ab, abpt); for (i = 0; i < sequences_.size(); ++i) free(bseqs[i]); free(bseqs); free(seq_lens);
        return true;
    } else { 
        abpoa_free(ab, abpt); for (i = 0; i < sequences_.size(); ++i) free(bseqs[i]); free(bseqs); free(seq_lens);
        return false;
    }
}

}
