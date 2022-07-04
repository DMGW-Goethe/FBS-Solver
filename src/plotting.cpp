#include "plotting.hpp"

#ifdef DEBUG_PLOTTING
namespace plt = matplotlibcpp;
#endif

void plotting::plot_evolution(const std::vector<integrator::step>& results, const std::vector<integrator::Event>& events, std::vector<int> plot_components, std::vector<std::string> labels, std::string filename) {
    assert(results.size() > 0);
    assert(plot_components.size() == labels.size());
#ifdef DEBUG_PLOTTING
    std::vector<double> r;
    r.reserve(results.size());
    for(auto it = results.begin(); it != results.end(); ++it)
        r.push_back(it->first);

    for(int i = 0; i < plot_components.size(); i++) {
        std::vector<double> y;
        y.reserve(results.size());
        int index = plot_components[i];
        assert(index >= 0 && index < results[0].second.size());
        for(auto it = results.begin(); it != results.end(); ++it)
            y.push_back(it->second[index]);
        plt::plot(r, y, {{"label", labels[i]}});
    }
    if(!filename.empty())
        plt::save(filename);
#else
    if(filename.empty())
        return;
    plotting::save_integration_data(results, plot_components, labels, filename.append(".txt"));
#endif
}

void plotting::save_integration_data(const std::vector<integrator::step>& results, std::vector<int> plot_components, std::vector<std::string> labels, std::string filename) {

	std::ofstream img;
    if(filename.empty())
        return;
	img.open(filename);

	if(img.is_open()) {
        img << "# r";
        for(int i = 0; i < plot_components.size(); i++)
            img << "\t" << labels[i];
		img << std::endl;

        for(auto it = results.begin(); it != results.end(); ++it) {
			img << std::fixed << std::setprecision(10) << it->first;    // radius
            for(int i = 0; i < plot_components.size(); i++)
                img << " " << it->second[plot_components[i]]; // the other variables
            img << std::endl;
		}
	}
	img.close();
}


void plotting::save_MR_data(std::vector<std::vector<double>> &array, std::vector<int> plot_components, std::vector<std::string> labels, std::string filename) {

    std::ofstream img;
    if(filename.empty())
        return;
	img.open(filename);

	if(img.is_open()) {
        // print the labels in line 1:
        for(int i = 0; i < plot_components.size(); i++)
            {img << labels[i] << "\t";}
		img << std::endl;

        // print all the data:
        for(unsigned i = 0; i < array.size(); ++i) {
            for(int j = 0; j < plot_components.size(); j++)
                {img << array[i][j] << " ";} // write the variables line by line
            img << std::endl;
		}
	}
	img.close();
}