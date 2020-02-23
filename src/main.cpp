/*


*/

#include <traversal.hpp>
#include <arg_parse.hpp>

#include <string>
#include <thread>
#include <vector>

#include <seqan3/range/views/char_to.hpp>
#include <seqan3/range/views/join.hpp>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/core/debug_stream.hpp>

#include <seqan3/argument_parser/all.hpp>

#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>

#include <range/v3/view/repeat.hpp>
#include <range/v3/core.hpp>
#include <range/v3/view/concat.hpp>
#include <range/v3/view/repeat_n.hpp>
#include <range/v3/all.hpp>

#include <seqan3/std/filesystem>

int main(int argc, char ** argv)
{
    auto args = initialize_argument_parser("HOPA", argc, argv);

    seqan3::sequence_file_input fin{args.in_file_path};


    auto config =   seqan3::align_cfg::mode{seqan3::global_alignment} | 
                    seqan3::align_cfg::aligned_ends{seqan3::free_ends_all} |
                    seqan3::align_cfg::scoring{seqan3::nucleotide_scoring_scheme
                        {
                            seqan3::match_score{4},
                            seqan3::mismatch_score{-2}
                        }
                    } |
                    seqan3::align_cfg::gap{seqan3::gap_scheme{seqan3::gap_score{-4}}} |
                    seqan3::align_cfg::parallel{args.threads} |
                    seqan3::align_cfg::result{seqan3::with_alignment};

// Ranges and views are the C++ equivalent of deboning a duck.

    // Get the top sequence: orient everything relative to this.            
    auto first_rec = *fin.begin();
    // If the user wants the orientation reversed, 
    // just reverse the first sequence at the beginning
    if(args.reverse)
        seqan3::get<seqan3::field::seq>(first_rec)=seqan3::views::complement(std::views::reverse(seqan3::get<seqan3::field::seq>(first_rec)));

    auto first_seq = seqan3::get<seqan3::field::seq>(first_rec);

    // Get all the other sequences in the file.
    auto back_recs = fin | seqan3::views::drop(1) | ranges::to<std::vector>();
    auto back_seqs = back_recs | std::views::transform([] (auto s) { return seqan3::get<seqan3::field::seq>(s) ;});
    // Duplicate each sequence
    auto back_seqs_2 = ranges::views::for_each(back_seqs, [] (auto c) {
        return ranges::yield_from(ranges::views::repeat_n(c,2));
    }) | ranges::to<std::vector>();
    // Couldn't figure out the range function, so got halfway and
    // went back to traditional iterators.
    for(auto it = std::begin(back_seqs_2); it < std::end(back_seqs_2); it++)
    {
        *(++it) = seqan3::views::complement(std::views::reverse(*it));
    }
    // And back to ranges!
    auto pair_seqs = back_seqs_2    | std::views::transform([first_seq](auto s){return std::pair{first_seq, s};})
                                    | ranges::to<std::vector>();

    auto results = seqan3::align_pairwise(pair_seqs, config);

    auto finit = std::begin(back_recs); 
    // Output the first sequence, which everything else is aligned back to.
    seqan3::sequence_file_output fout{args.out_file_path};
    auto temp = *finit;
    fout.push_back(first_rec);
    for (auto r = std::begin(results); r != std::end(results); r++)
    {
        auto forward_res = *r; 
        auto reverse_res = *(++r);
        if(forward_res.score() > reverse_res.score())
        {
            // If the forward sequence aligns better
            // then output the forward sequence information.
            fout.push_back(*finit);
        }
        else
        {
            // If the reverse sequence aligns better
            // then swap the forward sequence out and output.
            seqan3::get<seqan3::field::seq>(*finit)=seqan3::views::complement(std::views::reverse(seqan3::get<seqan3::field::seq>(*finit)));
            fout.push_back(*finit);
        }
        finit++;
    }
    return 0;
}
