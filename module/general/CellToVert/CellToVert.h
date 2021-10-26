#ifndef CELLTOVERT_H
#define CELLTOVERT_H

#include <vistle/module/module.h>
#include <vistle/core/vector.h>

class CellToVert: public vistle::Module {
public:
    CellToVert(const std::string &name, int moduleID, mpi::communicator comm);
    ~CellToVert();

private:
    static const int NumPorts = 3;

    std::vector<vistle::Port *> m_data_in, m_data_out;

    bool compute(std::shared_ptr<vistle::PortTask> task) const override;
};

#endif
