#ifndef NEWTONMETHOD_HPP
#define NEWTONMETHOD_HPP

#include <functional>
#include <utility>

class NewtonMethod
{
    public:

        NewtonMethod() {}

        virtual ~NewtonMethod() {}

        double solve(std::function<double(double)> f, std::function<double(double)> fd, double guess) const;

        std::pair<double, double> solve(std::function<double(double, double)> f, std::function<double(double, double)> fd1, std::function<double(double, double)> fd2,
                                        std::function<double(double, double)> g, std::function<double(double, double)> gd1, std::function<double(double, double)> gd2,
                                        double guess1, double guess2) const;

    private:
        static constexpr double numTol = 0.0001;

};

#endif // NEWTONMETHOD_HPP