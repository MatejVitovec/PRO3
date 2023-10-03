#ifndef CONTROL_HPP
#define CONTROL_HPP

#include <string>
#include <vector>
#include <map>

class Control
{
    public:

        Control() {}

        void readControl();
        void readBoundaryCondition();
        void parse();


        virtual ~Control() {}

    private:
        std::map<std::string, std::string> setupString;
        std::map<std::string, std::map<std::string, std::string>> boundaryConditionString;
        
};



#endif // CONTROL_HPP