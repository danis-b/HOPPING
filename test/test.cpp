#include <iostream>
#include "xtensor/xarray.hpp"
#include "xtensor/xio.hpp"
#include "xtensor/xview.hpp"

int main()
{
  xt::xarray<int> arr{1, 2, 3, 4, 5, 6, 7, 8, 9};
  

  arr1.reshape({3, 3});

  std::cout << arr1 << std::endl;
}