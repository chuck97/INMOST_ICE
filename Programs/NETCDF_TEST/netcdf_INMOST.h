#pragma once
#include "config.h"



class Dim
{
public:
    Dim(int id_,
        const std::string& name_,
        size_t len_);

    int GetId() const;
    std::string GetName() const;
    size_t GetSize() const;
    void DimInfo() const;

private:
    const int id;
    const std::string name;
    const size_t len;
};

class anyAtt
{
public:
    anyAtt(int id_,
           int owner_var_id_,
           const std::string& name_);

    virtual void AttInfo() const
    {};

protected:
    const int id;
    const int owner_var_id;
    const std::string name;
};

template<typename MessageType>
class Att: public anyAtt
{
public:
    Att(int id_,
        int owner_var_id_,
        const std::string& name_,
        const MessageType& message_):
        anyAtt(id_, owner_var_id_, name_),
        message(message_)
        { 
        };

    void AttInfo() const
    {
        std::cout << "Attribute '" << name
                  << "': id = " << id << ", "
                  << "owner var id = " << owner_var_id << ", "
                  << "message = " << message << std::endl;
    }

private:
    const MessageType message;
};

class anyVar
{
public:
    anyVar(const int& id_,
           const std::string& name_,
           const std::vector<Dim*>& dims_,
           const std::vector<std::shared_ptr<anyAtt>>& atts_);

    virtual void VarInfo() const
    {
    }

    virtual std::any Get_array1D() const
    {};

    virtual std::any Get_array2D() const
    {};

protected:
    const int id;
    const std::string name;
    const std::vector<Dim*> dims;
    const std::vector<std::shared_ptr<anyAtt>> atts;
};

template<typename VarType>
class Var : public anyVar
{
public: 
    Var(const int& _id,
        const std::string& _name,
        const std::vector<Dim*>& _dims,
        const std::vector<std::shared_ptr<anyAtt>>& _atts,
        nc_type type_,
        const std::vector<VarType>& _array1D, 
        const std::vector<std::vector<VarType>>& _array2D):
        anyVar(_id, _name, _dims, _atts),
        type(type_),
        array1D(_array1D),
        array2D(_array2D)
    {};

    void VarInfo() const
    {
        size_t lin_size;
        size_t x_size, y_size;
        if (array1D.empty())
        {
            lin_size = 0;
            if (array2D.empty())
            {
                x_size = 0;
                y_size = 0;
            }
            else
            {
                x_size = array2D.size();
                y_size = array2D.back().size();
            }
        }
        else
        {
            lin_size = array1D.size();
            x_size = 0;
            y_size = 0;
        }
        
        std::cout << "Variable: '" << name
                  << "': id = " << id << ", "
                  << "type = " << type  << ", "
                  << "ndims = " << dims.size()  << ", "
                  << "natts = " << atts.size()  << ", "
                  << "1d_size = "  << lin_size << ", "
                  << "2d_size = "  << x_size << "x" << y_size
                  << std::endl;
    }
    
    std::any Get_array1D() const
    {
        return array1D;
    }

    std::any Get_array2D() const
    {
        return array1D;
    }

private:
    const nc_type type;
    const std::vector<VarType> array1D;
    const std::vector<std::vector<VarType>> array2D;
};


class Dataset
{
public:
    Dataset(char* filename);
    void ReadDims();                          // +
    Dim* GetDimById(int id);                  // +
    void ReadAtts();                          // +
    void ReadVars(size_t time_step);          // + (Serial) - (Parallel)
    void Read();                              // -
    void DimsInfo() const;                    // +
    void AttsInfo() const;                    // +
    void VarsInfo() const;                    // +
    void FileInfo() const;                    // -
    void Write(const std::string& filename);  // -

private:
    int file_id;  
    std::vector<Dim> dims;
    std::vector<std::shared_ptr<anyAtt>> atts;
    std::vector<std::shared_ptr<anyVar>> vars;
    std::vector<std::shared_ptr<anyVar>> coords;
};