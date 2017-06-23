#include <iostream>
#include <string>
#include <fstream>
#include "kmer.h"
#include <vector>

bool readgenome(const std::string& file_name, std::string& genome) {
	std::ifstream in_str(file_name);
	if (!in_str.good()) {
		std::cerr << "Can't open " << file_name << " to read.\n";
		return false;
	}
	char tmp;
	while (in_str >> tmp) {
		genome.push_back(tmp);
	}
	return true;
}

unsigned int hash(std::string const& key) {        		
	unsigned int hash = 1315423911;
	for (unsigned int i = 0; i < key.length(); i++)
		hash ^= ((hash << 5) + key[i] + (hash >> 2));
	return hash;
}

bool create_table(const std::string& genome, std::vector<Kmer>& K_table, const unsigned int& kmer, const float& occupancy) {
	float current_occupancy = 0.0;
	unsigned int unique_keys = 0;
	unsigned int hashvalue_index;
	for (unsigned int i = 0; i + (kmer - 1) < genome.size(); ++i ) {
		hashvalue_index = hash(genome.substr(i, kmer)) % K_table.size();
		//Check if index has value
		if (K_table[hashvalue_index].get_Sequence() == "empty") {
			K_table[hashvalue_index].set_Sequence(genome.substr(i, kmer));
			K_table[hashvalue_index].insert_Postion(i);
			//Recalculate Occupany
			unique_keys += 1;
			current_occupancy = unique_keys / K_table.size();
			if (current_occupancy >= occupancy) {
				return false;
			}
		}
		//Kmer is already in table, add new postion
		else if (K_table[hashvalue_index].get_Sequence() == genome.substr(i, kmer)) {
			K_table[hashvalue_index].insert_Postion(i);
		}
		//Collsion - Will use Quadratic Probing
		else {
			for (int j = 1 ; j < K_table.size() ; ++j ) {
				//Back to orginal, table is full
				if (( hashvalue_index + j * j) % K_table.size() == hashvalue_index )
					return false;
				if (K_table[(hashvalue_index + j * j) % K_table.size()].get_Sequence() == genome.substr(i, kmer)) { 
					K_table[(hashvalue_index + j * j) % K_table.size()].insert_Postion(i);
					break;
				}
				else if ( K_table[(hashvalue_index + j * j) % K_table.size()].get_Sequence() == "empty") {
					K_table[(hashvalue_index + j * j) % K_table.size()].set_Sequence(genome.substr(i, kmer));
					K_table[(hashvalue_index + j * j) % K_table.size()].insert_Postion(i);
					//Recalculate Occupany
					unique_keys += 1;
					current_occupancy = unique_keys / K_table.size();
					if (current_occupancy >= occupancy) {
						return false;
					}
					break;
				}
			}
		}
	}
	return true;
}

void find(const std::string& genome, const std::vector<Kmer>& K_table, const unsigned int& kmer , 
	const unsigned int& mismatch, const std::string& query_string) {
	std::string query_kmer = query_string.substr(0, kmer);
	unsigned int hashvalue_index = hash(query_kmer) % K_table.size();
	unsigned int current_mismatch;
	std::string compare;
	bool first = false;
	if ( K_table[hashvalue_index].get_Sequence() == query_kmer) {
		// Same Kmer has shown up mutliple times
		if (K_table[hashvalue_index].get_Postions().size() >= 2) {
			for (int i = 0; i < K_table[hashvalue_index].get_Postions().size(); ++i) {
				current_mismatch = 0;
				//Check to stay in bounds
				if (K_table[hashvalue_index].get_Postions()[i] + (query_string.size() - 1) >= genome.size() ) {
					continue;
				}
				// Compare Values
				compare = genome.substr(K_table[hashvalue_index].get_Postions()[i] , query_string.size() );
				for ( int j = 0; j < compare.size() ; ++j) {
					if (compare[j] != query_string[j]) {
						current_mismatch += 1;
						if (current_mismatch > mismatch)
							break;
					}
				}
				//Found Word
				if (current_mismatch <= mismatch) {
					if (!first ) {
						std::cout << "Query: " << query_string << '\n' << K_table[hashvalue_index].get_Postions()[i] << 
						" " << current_mismatch << " " << compare << std::endl;
						first = true;
					}
					else
						std::cout << K_table[hashvalue_index].get_Postions()[i] << " " << current_mismatch << " " << compare << std::endl;
				}
			}
			// Not found
			if (!first)
				std::cout << "Query: " << query_string << '\n' << "No Match" << std::endl;
			return;
		}
		//Make sure we stay in bounds for the rest of the operations
		if (K_table[hashvalue_index].get_Postions()[0] + (query_string.size() - 1) >= genome.size()) {
			std::cout << "Query: " << query_string << '\n' << "No Match" << std::endl;
			return;
		}
		//Single Occurance of Kmer
		compare = genome.substr(K_table[hashvalue_index].get_Postions()[0] , query_string.size() );
		current_mismatch = 0;
		for ( int j = 0; j < compare.size() ; ++j) {
			if (compare[j] != query_string[j]) {
				current_mismatch += 1;
				if (current_mismatch > mismatch) {
					std::cout << "Query: " << query_string << '\n' << "No Match" << std::endl;
					return;
				}
			}
		}
		std::cout << "Query: " << query_string << '\n' << K_table[hashvalue_index].get_Postions()[0] <<
		 " " << current_mismatch << " " << compare << std::endl;
		return;
	}
	//Check for Collision
	for (int i = 1 ; i < K_table.size() ; ++i ) {
		// Table is full
		if ((hashvalue_index + i * i) % K_table.size() == hashvalue_index) {
			std::cout << "Query: " << query_string << '\n' << "No Match" << std::endl;
			return;
		}
		if (K_table[(hashvalue_index + i * i) % K_table.size()].get_Sequence() == query_kmer) {
			// Same Kmer has shown up mutliple times
			if (K_table[(hashvalue_index + i * i) % K_table.size()].get_Postions().size() >= 2) {
				for (int j = 0; j < K_table[(hashvalue_index + i * i) % K_table.size()].get_Postions().size(); ++j) {
					current_mismatch = 0;
					//Check to stay in bounds
					if (K_table[(hashvalue_index + i * i) % K_table.size()].get_Postions()[j] + (query_string.size() - 1) < genome.size() ) {
						compare = genome.substr(K_table[(hashvalue_index + i * i) % K_table.size()].get_Postions()[j] , query_string.size() );
						for ( int n = 0; n < compare.size() ; ++n) {
							if (compare[n] != query_string[n]) {
								current_mismatch += 1;
								if (current_mismatch > mismatch)
									break;
							}
						}
						if (current_mismatch > mismatch)
							continue;
						// Found
						if (!first) {
							std::cout << "Query: " << query_string << '\n' << K_table[(hashvalue_index + i * i) % K_table.size()].get_Postions()[j] <<
							 " " << current_mismatch << " " << compare << std::endl;
							first = true;
						}
						else
							std::cout << K_table[(hashvalue_index + i * i) % K_table.size()].get_Postions()[j] << " " << current_mismatch <<
							 " " << compare << std::endl;
					}
				}
				//Not Found
				if (!first)
					std::cout << "Query: " << query_string << '\n' << "No Match" << std::endl;
				return;
			}
			//Single Occurance of Kmer && make sure we stay in bounds
			if (K_table[(hashvalue_index + i * i) % K_table.size()].get_Postions()[0] + (query_string.size() - 1) >= genome.size()) {
				std::cout << "Query: " << query_string << '\n' << "No Match" << std::endl;
				return;
			}
			compare = genome.substr(K_table[(hashvalue_index + i * i) % K_table.size()].get_Postions()[0] , query_string.size() );
			current_mismatch = 0;
			for ( int j = 0; j < compare.size() ; ++j) {
				if (compare[j] != query_string[j]) {
					current_mismatch += 1;
					if (current_mismatch > mismatch) {
						std::cout << "Query: " << query_string << '\n' << "No Match" << std::endl;
						return;
					}
				}
			}
			std::cout << "Query: " << query_string << '\n' << K_table[(hashvalue_index + i * i) % K_table.size()].get_Postions()[0] << 
			" " << current_mismatch << " " << compare << std::endl;
			return;
		}
	}
	//Not found
	std::cout << "Query: " << query_string << '\n' << "No Match" << std::endl;
}

int main() {
	std::string genome; std::string file_name; std::string query_string; std::string tmp;
	unsigned int kmer; unsigned int table_size = 100; unsigned int mismatch;
	float occupancy = 0.5;
	bool constructed = false;
	std::vector<Kmer> K_table;
	//Start processes
	std::cin >> tmp;
	while ( tmp != "quit") {
		if (tmp == "genome") {
			std::cin >> file_name;
			if (!readgenome(file_name, genome))
				return 1;
		}
		if (tmp == "table_size") {
			std::cin >> table_size;
			//Set length of table
			if (!constructed) {
				K_table.resize(table_size);
				constructed = true;
			}
		}
		if (tmp == "occupancy") {
			std::cin >> occupancy;
		}
		if (tmp == "kmer") {
			std::cin >> kmer;
			if (!constructed) {
				K_table.resize(table_size);
				constructed = true;
			}
			//Check if we need to resize the table, occupancy failed
			while (!create_table(genome, K_table, kmer, occupancy)) {
				table_size = 2 * table_size;
				K_table.clear();
				K_table.resize(table_size);
			}
		}
		if (tmp == "query") {
			mismatch = 0;
			std::cin >> mismatch >> query_string;
			find(genome, K_table, kmer, mismatch, query_string);
		}
		std::cin >> tmp;
	}
	return 0;
}