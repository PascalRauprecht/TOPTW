#include <boost/lambda/lambda.hpp>
#include <vector>
#include <algorithm>
#include <iostream>
#include "instance.h"

int main()
{
    std::vector<int> v{1, 3, 2};
    std::for_each(v.begin(), v.end(),
                  std::cout << boost::lambda::_1 << "\n");
}

//test
//test2

//test7

//test11
