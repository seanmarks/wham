#include "DataSummary.h"


DataSummary::DataSummary()
{}


DataSummary::DataSummary(
	const std::string& data_summary_file,
	const ParameterPack& input_pack
):
	data_summary_file_(data_summary_file),
	data_summary_path_( FileSystem::get_basename(data_summary_file_) ),  // full path
	col_data_label_(0),
	col_t_min_(3), col_t_max_(col_t_min_+1), col_T_(-1)
{
	// Largest column index of interest
	int max_col = std::max( { col_data_label_, col_t_min_, col_T_ }, std::less<int>() );

	std::ifstream ifs(data_summary_file_);
	if ( not ifs.is_open() ) {
		throw std::runtime_error("Unable to open data summary file: " + data_summary_file_ + "\n");
	}

	std::string line, token;
	std::vector<std::string> tokens;

	while ( getline(ifs, line) ) {
		std::stringstream ss(line);
		ss >> token;

		if ( line.empty() or token[0] == '#' ) {
			continue;  // skip empty lines and comments
		}
		
		// Tokenize the line
		tokens = {{ token }};
		while ( ss >> token ) {
			tokens.push_back( token );
		}

		// Check that the appropriate number of tokens was parsed
		int num_tokens = tokens.size();
		if ( max_col >= num_tokens ) {
			throw std::runtime_error("a line in the data summary has too few columns");
		}

		// Data set label
		data_set_labels_.push_back( tokens[col_data_label_] );

		// Production phase
		t_min_.push_back( std::stod(tokens[col_t_min_]) );
		t_max_.push_back( std::stod(tokens[col_t_max_]));

		/*
		// Temperature (optional)
		if ( col_T_ >= 0 ) {
			temperatures_.push_back( std::stod(tokens[col_T_]) );
		}
		*/
	}

	int num_simulations = data_set_labels_.size();
	for ( int i=0; i<num_simulations; ++i ) {
		auto ret = map_data_set_labels_to_indices_.insert( std::make_pair( data_set_labels_[i], i ) );
		if ( ret.second == false ) {
			throw std::runtime_error("Data set label \"" + data_set_labels_[i] + "\" was defined more than once");
		}
	}
}
