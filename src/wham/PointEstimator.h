#ifndef POINT_ESTIMATOR_H
#define POINT_ESTIMATOR_H

#include <vector>

#include "Statistics.h"

template<typename T>
class PointEstimator
{
 public:
	PointEstimator() {}

	PointEstimator(const int num_samples):
		samples_(0)
	{
		samples_.reserve(num_samples);
	}

	void clear() {
		samples_.clear();
	}

	void addSample(const T& sample) {
		samples_.push_back(sample);
	}

	// Statistics
	template<typename F>
	F average() const {
		return Statistics::average(samples_);
	}
	template<typename F>
	F std_dev(const int delta_dof = 1) const {
		return Statistics::std_dev(samples_, delta_dof);
	}

 private:
	std::vector<T> samples_;
};

#endif // ifndef POINT_ESTIMATOR_H
