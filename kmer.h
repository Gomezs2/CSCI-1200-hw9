#ifndef _kmer_h_
#define _kmer_h_

#include <vector>

class Kmer {
public:
	Kmer() : sequence_("empty"), postions_() {}
	Kmer(const std::string& sequence, unsigned int postion) { sequence_ = sequence; postions_.push_back(postion); }
	const std::string& get_Sequence() const {return sequence_;}
	const std::vector<unsigned int>& get_Postions() const {return postions_;}

	void set_Sequence(const std::string& sequence) {sequence_ = sequence;}
	void insert_Postion(unsigned int location) {postions_.push_back(location);}
private:
	std::string sequence_;
	std::vector<unsigned int> postions_;
};

#endif