# Histogram
Histogram in C++11, 1D, simple header-only, inspired by R, calculate and optimize breaks automatically. Accepts different precissions.
# Documentation

<!-- [Doxygen generated Docs](http://phcerdan.github.io/histo-header/docs/html/index.html) -->
[Doxygen generated Docs](http://phcerdan.github.io/histogram/docs/html/index.html)

Warning, docs are not generated automatically and might be outdated.

# Downloading
You can use the header as is, just download it, and include it in your project.
MIT License, but please contribute back if you add any extra feature!

# Using it.
You should read the test file as a sample for usage.
```cpp
vector<double> data{1.0, 1.0, 2.0, 3.0, 19.0, 0.5, 14.0};
histo::Histo<double> h_scott(data); // Default method to calculate Breaks.
size_t bins = 10;
auto breaks_with_bins = histo::GenerateBreaksFromRangeAndBins<double>(0.0, 20.0, bins);
double width = 1.0;
auto breaks_with_fixed_width = histo::GenerateBreaksFromRangeAndWidth<double>(0.0, 20.0, width);
histo::Histo<double> h_with_bins(data, breaks_with_bins);
histo::Histo<double> h_with_width(data, breaks_with_width);
```

Each histogram has public members: `bins`, `breaks`, `counts` and `range`.
We can fill the histogram with `FillCounts(data)`, called at constructor.
The data is not stored in the histogram.

We can fill the bins with more data to an existing histogram.
```cpp
vector<double> extra_data{7.0, 13.0};
h_with_bins.FillCounts(extra_data);
```

We can also normalize the histogram to get a probability density function from it.

```cpp
histo::Histo<double, unsigned int> regular_histo(data);
histo::Histo<double, double> normalized_histogram = histo::NormalizeByArea(regular_histo);
```

Optionally, we can use VTK (vtkChartXY) to visualize the histogram.

```cpp
vector<double> data{0.0, 1.0, 1.0,1.0, 2.0, 3.0, 5.0, 5.0, 8.0, 8.0,  12.0};
Histo<double> h(data, histo::GenerateBreaksFromRangeAndBins<double>(0.0,15.0, 5));
h.PrintBreaksAndCounts(std::cout);
h.name = "withJustData";
visualize_histo(h, vtkChart::LINE);
visualize_histo(h, vtkChart::BAR);
```


<img src="https://github.com/phcerdan/histogram/blob/gh-pages/readme_images/just_data_line.png" alt="Line" width="640" height="480">
<img src="https://github.com/phcerdan/histogram/blob/gh-pages/readme_images/bar_just_data.png" alt="Bar" width="640" height="480">

# Test
All the features are tested using gtest.
