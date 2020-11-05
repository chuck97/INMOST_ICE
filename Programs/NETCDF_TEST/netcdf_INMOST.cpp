#include "netcdf_INMOST.h"

Dim::Dim(int id_,
        const std::string& name_,
        size_t len_):
        id(id_),
        name(name_),
        len(len_)
{
}

int Dim::GetId() const
{
    return id;
}

std::string Dim::GetName() const
{
    return name;
}

size_t Dim::GetSize() const
{
    return len;
}

void Dim::DimInfo() const
{
    std::cout << "Dimension '" << name 
              << "': id = " << id << ", "
              << "size = " << len << std::endl;
}

anyAtt::anyAtt(int id_,
               int owner_var_id_,
               const std::string& name_):
               id(id_),
               owner_var_id(owner_var_id_),
               name(name_)
{
}


anyVar::anyVar(const int& id_,
               const std::string& name_,
               const std::vector<Dim*>& dims_,
               const std::vector<std::shared_ptr<anyAtt>>& atts_):
               id(id_),
               name(name_),
               dims(dims_),
               atts(atts_)
{
};


void Dataset::ReadDims()
{
    if (!dims.empty())
    {
        throw std::invalid_argument("Can't read dims twice");
        return;
    }
    int include_parents, ndims;
    int dimids[NC_MAX_DIMS+1];
    int retval;
    if ((retval = nc_inq_dimids(file_id, &ndims, dimids, include_parents)))
        ERR(retval);
    for (size_t i = 0; i < ndims; ++i)
    {
        char name[NC_MAX_NAME+1];
        size_t dimlen;
        if ((retval = nc_inq_dim(file_id, dimids[i], name, &dimlen)))
            ERR(retval);
        std::string nm = name;
        dims.push_back({(int)i, nm, dimlen});
    }
    std::cout << "Read Dimensions: successful" << std::endl;
}

void Dataset::ReadAtts()
{
    if (!atts.empty())
    {
        throw std::invalid_argument("Can't read atts twice");
        return;
    }

    int natts;
    int retval;

    if((retval = nc_inq(file_id, NULL, NULL, &natts, NULL)))
        ERR(retval);
    
    for (size_t i = 0; i < natts; ++i)
    {
        char name[NC_MAX_NAME+1];
        if((retval = nc_inq_attname(file_id, NC_GLOBAL, i, name)))
            ERR(retval)

        nc_type type;
        size_t len;

        if((retval = nc_inq_att(file_id, NC_GLOBAL, name, &type, &len)))
            ERR(retval)
        
        if (type == NC_CHAR)
        {
            char attvaluec[len+1];

            if((retval = nc_get_att(file_id, NC_GLOBAL, name, attvaluec)))
                ERR(retval)

            std::string attvalue(attvaluec);
            std::string attname(name);

            std::shared_ptr<anyAtt> stringatt;
            stringatt = std::make_shared<Att<std::string>>((int)i, NC_GLOBAL, attname, attvalue);
            atts.push_back(stringatt);
        }
        else if((type == NC_INT) or (type == NC_LONG) or (type == NC_SHORT) or (type == NC_BYTE))
        {
            int attvaluei;

            if((retval = nc_get_att(file_id, NC_GLOBAL, name, &attvaluei)))
                ERR(retval)

            std::string attname(name);

            std::shared_ptr<anyAtt> intatt;
            intatt = std::make_shared<Att<int>>((int)i, NC_GLOBAL, attname, attvaluei);
            atts.push_back(intatt);
        }
        else if(type == NC_DOUBLE)
        {
            double attvalued;

            if((retval = nc_get_att(file_id, NC_GLOBAL, name, &attvalued)))
                ERR(retval)

            std::string attname(name);

            std::shared_ptr<anyAtt> doubleatt;
            doubleatt = std::make_shared<Att<int>>((int)i, NC_GLOBAL, attname, attvalued);
            atts.push_back(doubleatt);
        }
        else if(type == NC_FLOAT)
        {
            float attvaluef;

            if((retval = nc_get_att(file_id, NC_GLOBAL, name, &attvaluef)))
                ERR(retval)

            std::string attname(name);

            std::shared_ptr<anyAtt> floatatt;
            floatatt = std::make_shared<Att<int>>((int)i, NC_GLOBAL, attname, attvaluef);
            atts.push_back(floatatt);
        }
        else if(type == NC_UINT)
        {
            uint attvalueu;

            if((retval = nc_get_att(file_id, NC_GLOBAL, name, &attvalueu)))
                ERR(retval)

            std::string attname(name);

            std::shared_ptr<anyAtt> uintatt;
            uintatt = std::make_shared<Att<int>>((int)i, NC_GLOBAL, attname, attvalueu);
            atts.push_back(uintatt);
        }
        else
        {
            throw std::invalid_argument("unknown attribute type: "+ std::to_string(type));
        }


    }
    std::cout << "Read Global Attributes: successful" << std::endl;
}

void Dataset::DimsInfo() const
{
    int i = 0;
    std::cout << "There are " << dims.size() << " dimensions:" << std::endl;
    for (const Dim& item : dims)
    {
        item.DimInfo();
    }
}

void Dataset::AttsInfo() const
{
    int i = 0;
    std::cout << "There are " << atts.size() << " global attributes:" << std::endl;
    for (const std::shared_ptr<anyAtt>& item : atts)
    {
        (*item).AttInfo();
    }
}



Dataset::Dataset(char* filename)
{
    int retval;
    int ncid;
    if ((retval = nc_open(filename, NC_NOWRITE, &ncid)))
        ERR(retval);
    file_id = ncid;
}

Dim* Dataset::GetDimById(int id)
{
    for (Dim& item: dims)
    {
        if (item.GetId() == id)
        {
            return &item;
        }
    }
}

void Dataset::ReadVars(size_t time_step)
{
    if (!vars.empty())
    {
        throw std::invalid_argument("Can't read vars twice");
        return;
    }

    if (dims.empty())
    {
        ReadDims();
    }    


    int nvars;
    int retval;

    if((retval = nc_inq(file_id, NULL, &nvars, NULL, NULL)))
        ERR(retval);
    
    for (size_t i = 0; i < nvars; ++i)
    {
        // Read var name
        char varnamec[NC_MAX_NAME+1];
        if((retval = nc_inq_var(file_id, i, varnamec, NULL, NULL, NULL, NULL)))
            ERR(retval);

        std::string varname(varnamec);

        // Read var type
        nc_type type;
        if((retval = nc_inq_var(file_id, i, NULL, &type, NULL, NULL, NULL)))
            ERR(retval);
    
        // Read var dims
        int nvardims;
        int vardimids[NC_MAX_DIMS + 1];
        bool time_dep = false;

        if((retval = nc_inq_var(file_id, i, NULL, NULL, &nvardims, vardimids, NULL)))
            ERR(retval);
        
        std::vector<Dim*> current_var_dims;
        for (size_t j = 0; j < nvardims; ++j)
        {
            Dim* d = GetDimById(vardimids[j]);
            if (d->GetName() == TIME_NAME)
            {
                time_dep = true;
            }
            current_var_dims.push_back(GetDimById(vardimids[j]));
        }
        
        // Read var attributes
        std::vector<std::shared_ptr<anyAtt>> current_var_atts;
        int natts;
        if((retval = nc_inq_var(file_id, (int)i, NULL, NULL, NULL, NULL, &natts)))
            ERR(retval);
        
        for (size_t j = 0; j < natts; ++j)
        {
            
            char name[NC_MAX_NAME+1];
            if((retval = nc_inq_attname(file_id, (int)i, (int)j, name)))
                ERR(retval)
            nc_type type;
            size_t len;


            if((retval = nc_inq_att(file_id, (int)i, name, &type, &len)))
                ERR(retval)

            if (type == NC_CHAR)
            {
                char attvaluec[len+1];

                if((retval = nc_get_att(file_id, (int)i, name, attvaluec)))
                    ERR(retval)

                std::string attvalue(attvaluec);
                std::string attname(name);

                std::shared_ptr<anyAtt> stringatt;
                stringatt = std::make_shared<Att<std::string>>((int)j, i, attname, attvalue);
                current_var_atts.push_back(stringatt);
            }
            else if((type == NC_INT) or (type == NC_LONG) or (type == NC_SHORT) or (type == NC_BYTE))
            {
                int attvaluei;

                if((retval = nc_get_att(file_id, (int)i, name, &attvaluei)))
                    ERR(retval)

                std::string attname(name);

                std::shared_ptr<anyAtt> intatt;
                intatt = std::make_shared<Att<int>>((int)j, (int)i, attname, attvaluei);
                current_var_atts.push_back(intatt);
            }
            else if(type == NC_DOUBLE)
            {
                double attvalued;

                if((retval = nc_get_att(file_id, (int)i, name, &attvalued)))
                    ERR(retval)

                std::string attname(name);

                std::shared_ptr<anyAtt> doubleatt;
                doubleatt = std::make_shared<Att<int>>((int)j, (int)i, attname, attvalued);
                current_var_atts.push_back(doubleatt);
            }
            else if(type == NC_FLOAT)
            {
                float attvaluef;

                if((retval = nc_get_att(file_id, (int)i, name, &attvaluef)))
                    ERR(retval)

                std::string attname(name);

                std::shared_ptr<anyAtt> floatatt;
                floatatt = std::make_shared<Att<int>>((int)j, (int)i, attname, attvaluef);
                current_var_atts.push_back(floatatt);
            }
            else if(type == NC_UINT)
            {
                uint attvalueu;

                if((retval = nc_get_att(file_id, (int)i, name, &attvalueu)))
                    ERR(retval)

                std::string attname(name);

                std::shared_ptr<anyAtt> uintatt;
                uintatt = std::make_shared<Att<int>>((int)j, (int)i, attname, attvalueu);
                current_var_atts.push_back(uintatt);
            }
            else
            {
                throw std::invalid_argument("unknown attribute type: "+ std::to_string(type));
            }
        }
         
        // Read Vars values
        
        if ((nvardims == 1) and !time_dep) // 1-d values
        {
            size_t varsize = current_var_dims.back()->GetSize();
            nc_type vartype;
            if ((retval = nc_inq_var(file_id, (int)i, NULL, &vartype, NULL, NULL, NULL)))
                ERR(retval)

            if (vartype == NC_FLOAT)
            {
                float values[varsize];
                if ((retval = nc_get_var_float(file_id, (int)i, values)))
                    ERR(retval)
                std::vector<float> float_array1D;
                for (size_t k = 0; k < varsize; ++k)
                {
                    float_array1D.push_back(values[k]);
                }
                std::vector<std::vector<float>> float_array2D{{}};

                std::shared_ptr<anyVar> floatvar;
                floatvar = std::make_shared<Var<float>>((int)i, varname, current_var_dims,
                                             current_var_atts, vartype, 
                                             float_array1D,
                                             float_array2D);
                vars.push_back(floatvar);
            }
            else if (vartype == NC_DOUBLE)
            {
                double values[varsize];
                if ((retval = nc_get_var_double(file_id, (int)i, values)))
                    ERR(retval)
                std::vector<double> double_array1D;
                for (size_t k = 0; k < varsize; ++k)
                {
                    double_array1D.push_back(values[k]);
                }
                std::vector<std::vector<double>> double_array2D{{}};
                std::shared_ptr<anyVar> doublevar;
                doublevar = std::make_shared<Var<double>>((int)i, varname, current_var_dims,
                                             current_var_atts, vartype, 
                                             double_array1D,
                                             double_array2D);
                vars.push_back(doublevar);
            }
            else if ((vartype == NC_INT) or (vartype == NC_LONG) or (vartype == NC_SHORT) or (vartype == NC_BYTE))
            {
                int values[varsize];
                if ((retval = nc_get_var(file_id, (int)i, values)))
                    ERR(retval)
                std::vector<int> int_array1D;
                for (size_t k = 0; k < varsize; ++k)
                {
                    int_array1D.push_back(values[k]);
                }
                std::vector<std::vector<int>> int_array2D{{}};
                std::shared_ptr<anyVar> intvar;
                intvar = std::make_shared<Var<int>>((int)i, varname, current_var_dims,
                                             current_var_atts, vartype, 
                                             int_array1D,
                                             int_array2D);
                vars.push_back(intvar);
            }
            else if (vartype == NC_UINT)
            {
                uint values[varsize];
                if ((retval = nc_get_var(file_id, (int)i, values)))
                    ERR(retval)
                std::vector<uint> uint_array1D;
                for (size_t k = 0; k < varsize; ++k)
                {
                    uint_array1D.push_back(values[k]);
                }
                std::vector<std::vector<uint>> uint_array2D{{}};
                std::shared_ptr<anyVar> uintvar;
                uintvar = std::make_shared<Var<uint>>((int)i, varname, current_var_dims,
                                             current_var_atts, vartype, 
                                             uint_array1D,
                                             uint_array2D);
                vars.push_back(uintvar);
            }
            else
            {
                throw std::invalid_argument("unknown attribute type: "+ std::to_string(vartype));
            }
        }
        else if((nvardims == 2) and !time_dep) // 2-d values without time
        {
            size_t size_x, size_y;
            if (current_var_dims[0]->GetName() == X_NAME)
            {
                size_x = current_var_dims[0]->GetSize();
                size_y = current_var_dims[1]->GetSize();
            }
            else
            {
                size_x = current_var_dims[1]->GetSize();
                size_y = current_var_dims[0]->GetSize();
            }
            
            nc_type vartype;
            if ((retval = nc_inq_var(file_id, (int)i, NULL, &vartype, NULL, NULL, NULL)))
                ERR(retval)

            if (vartype == NC_FLOAT)
            {
                //////////////////////////////////////////////////////////////////////
                //////////////////////////////////////////////////////////////////////
                float values[size_y][size_x];
                static size_t start[] = {0, 0};
                static size_t count[] = {size_y, size_x};
                if ((retval = nc_get_vara_float(file_id, (int)i, start, count, &values[0][0])))
                    ERR(retval)
                std::vector<std::vector<float>> float_array2D;
                for (size_t k_y = 0; k_y < size_y; ++k_y)
                {
                    std::vector<float> float_array_row;
                    for (size_t k_x = 0; k_x < size_x; ++k_x)
                    {
                        float_array_row.push_back(values[k_y][k_x]);
                    }
                    float_array2D.push_back(float_array_row);
                }
                //////////////////////////////////////////////////////////////////////
                //////////////////////////////////////////////////////////////////////
                std::vector<float> float_array1D{};

                std::shared_ptr<anyVar> floatvar;
                floatvar = std::make_shared<Var<float>>((int)i, varname, current_var_dims,
                                             current_var_atts, vartype, 
                                             float_array1D,
                                             float_array2D);
                vars.push_back(floatvar);
            }
            else if (vartype == NC_DOUBLE)
            {
                //////////////////////////////////////////////////////////////////////
                //////////////////////////////////////////////////////////////////////
                double values[size_y][size_x];
                static size_t start[] = {0, 0};
                static size_t count[] = {size_y, size_x};
                if ((retval = nc_get_vara_double(file_id, (int)i, start, count, &values[0][0])))
                    ERR(retval)
                std::vector<std::vector<double>> double_array2D;
                for (size_t k_y = 0; k_y < size_y; ++k_y)
                {
                    std::vector<double> double_array_row;
                    for (size_t k_x = 0; k_x < size_x; ++k_x)
                    {
                        double_array_row.push_back(values[k_y][k_x]);
                    }
                    double_array2D.push_back(double_array_row);
                }
                //////////////////////////////////////////////////////////////////////
                //////////////////////////////////////////////////////////////////////
                std::vector<double> double_array1D{};

                std::shared_ptr<anyVar> doublevar;
                doublevar = std::make_shared<Var<double>>((int)i, varname, current_var_dims,
                                             current_var_atts, vartype, 
                                             double_array1D,
                                             double_array2D);
                vars.push_back(doublevar);
            }
            else if ((vartype == NC_INT) or (vartype == NC_LONG) or (vartype == NC_SHORT) or (vartype == NC_BYTE))
            {
                //////////////////////////////////////////////////////////////////////
                //////////////////////////////////////////////////////////////////////
                int values[size_y][size_x];
                static size_t start[] = {0, 0};
                static size_t count[] = {size_y, size_x};
                if ((retval = nc_get_vara_int(file_id, (int)i, start, count, &values[0][0])))
                    ERR(retval)
                std::vector<std::vector<int>> int_array2D;
                for (size_t k_y = 0; k_y < size_y; ++k_y)
                {
                    std::vector<int> int_array_row;
                    for (size_t k_x = 0; k_x < size_x; ++k_x)
                    {
                        int_array_row.push_back(values[k_y][k_x]);
                    }
                    int_array2D.push_back(int_array_row);
                }
                //////////////////////////////////////////////////////////////////////
                //////////////////////////////////////////////////////////////////////
                std::vector<int> int_array1D{};

                std::shared_ptr<anyVar> intvar;
                intvar = std::make_shared<Var<int>>((int)i, varname, current_var_dims,
                                             current_var_atts, vartype, 
                                             int_array1D,
                                             int_array2D);
                vars.push_back(intvar);
            }
            else if (vartype == NC_UINT)
            {
                //////////////////////////////////////////////////////////////////////
                //////////////////////////////////////////////////////////////////////
                uint values[size_y][size_x];
                static size_t start[] = {0, 0};
                static size_t count[] = {size_y, size_x};
                if ((retval = nc_get_vara_uint(file_id, (int)i, start, count, &values[0][0])))
                    ERR(retval)
                std::vector<std::vector<uint>> uint_array2D;
                for (size_t k_y = 0; k_y < size_y; ++k_y)
                {
                    std::vector<uint> uint_array_row;
                    for (size_t k_x = 0; k_x < size_x; ++k_x)
                    {
                        uint_array_row.push_back(values[k_y][k_x]);
                    }
                    uint_array2D.push_back(uint_array_row);
                }
                //////////////////////////////////////////////////////////////////////
                //////////////////////////////////////////////////////////////////////
                std::vector<uint> uint_array1D{};

                std::shared_ptr<anyVar> uintvar;
                uintvar = std::make_shared<Var<uint>>((int)i, varname, current_var_dims,
                                             current_var_atts, vartype, 
                                             uint_array1D,
                                             uint_array2D);
                vars.push_back(uintvar);
            }
            else
            {
                throw std::invalid_argument("unknown variable type: "+ std::to_string(vartype));
            }
        }
        else if((nvardims == 3) and time_dep) // 3-d values with time
        {
            size_t size_x, size_y;
            if (current_var_dims[0]->GetName() == TIME_NAME)
            {
                if (current_var_dims[1]->GetName() == X_NAME)
                {
                    size_x = current_var_dims[1]->GetSize();
                    size_y = current_var_dims[2]->GetSize();
                }
                else
                {
                    size_x = current_var_dims[2]->GetSize();
                    size_y = current_var_dims[1]->GetSize();
                }
            }
            else if (current_var_dims[1]->GetName() == TIME_NAME)
            {
                if (current_var_dims[0]->GetName() == X_NAME)
                {
                    size_x = current_var_dims[0]->GetSize();
                    size_y = current_var_dims[2]->GetSize();
                }
                else
                {
                    size_x = current_var_dims[2]->GetSize();
                    size_y = current_var_dims[0]->GetSize();
                }
            }
            else
            {
                if (current_var_dims[0]->GetName() == X_NAME)
                {
                    size_x = current_var_dims[0]->GetSize();
                    size_y = current_var_dims[1]->GetSize();
                }
                else
                {
                    size_x = current_var_dims[1]->GetSize();
                    size_y = current_var_dims[0]->GetSize();
                }
            }
            
            nc_type vartype;
            if ((retval = nc_inq_var(file_id, (int)i, NULL, &vartype, NULL, NULL, NULL)))
                ERR(retval)
            
            if (vartype == NC_FLOAT)
            {
                //////////////////////////////////////////////////////////////////////
                //////////////////////////////////////////////////////////////////////
                float values[size_y][size_x];
                static size_t start[] = {time_step, 0, 0};
                static size_t count[] = {1, size_y, size_x};
                if ((retval = nc_get_vara_float(file_id, (int)i, start, count, &values[0][0])))
                    ERR(retval)
                std::vector<std::vector<float>> float_array2D;
                for (size_t k_y = 0; k_y < size_y; ++k_y)
                {
                    std::vector<float> float_array_row;
                    for (size_t k_x = 0; k_x < size_x; ++k_x)
                    {
                        float_array_row.push_back(values[k_y][k_x]);
                    }
                    float_array2D.push_back(float_array_row);
                }
                //////////////////////////////////////////////////////////////////////
                //////////////////////////////////////////////////////////////////////
                std::vector<float> float_array1D{};

                std::shared_ptr<anyVar> floatvar;
                floatvar = std::make_shared<Var<float>>((int)i, varname, current_var_dims,
                                             current_var_atts, vartype, 
                                             float_array1D,
                                             float_array2D);
                vars.push_back(floatvar);
            }
            else if (vartype == NC_DOUBLE)
            {
                //////////////////////////////////////////////////////////////////////
                //////////////////////////////////////////////////////////////////////
                double values[size_y][size_x];
                static size_t start[] = {time_step, 0, 0};
                static size_t count[] = {1, size_y, size_x};
                if ((retval = nc_get_vara_double(file_id, (int)i, start, count, &values[0][0])))
                    ERR(retval)
                std::vector<std::vector<double>> double_array2D;
                for (size_t k_y = 0; k_y < size_y; ++k_y)
                {
                    std::vector<double> double_array_row;
                    for (size_t k_x = 0; k_x < size_x; ++k_x)
                    {
                        double_array_row.push_back(values[k_y][k_x]);
                    }
                    double_array2D.push_back(double_array_row);
                }
                //////////////////////////////////////////////////////////////////////
                //////////////////////////////////////////////////////////////////////
                std::vector<double> double_array1D{};

                std::shared_ptr<anyVar> doublevar;
                doublevar = std::make_shared<Var<double>>((int)i, varname, current_var_dims,
                                             current_var_atts, vartype, 
                                             double_array1D,
                                             double_array2D);
                vars.push_back(doublevar);
            }
            else if ((vartype == NC_INT) or (vartype == NC_LONG) or (vartype == NC_SHORT) or (vartype == NC_BYTE))
            {
                //////////////////////////////////////////////////////////////////////
                //////////////////////////////////////////////////////////////////////
                int values[size_y][size_x];
                static size_t start[] = {time_step, 0, 0};
                static size_t count[] = {1, size_y, size_x};
                if ((retval = nc_get_vara_int(file_id, (int)i, start, count, &values[0][0])))
                    ERR(retval)
                std::vector<std::vector<int>> int_array2D;
                for (size_t k_y = 0; k_y < size_y; ++k_y)
                {
                    std::vector<int> int_array_row;
                    for (size_t k_x = 0; k_x < size_x; ++k_x)
                    {
                        int_array_row.push_back(values[k_y][k_x]);
                    }
                    int_array2D.push_back(int_array_row);
                }
                //////////////////////////////////////////////////////////////////////
                //////////////////////////////////////////////////////////////////////
                std::vector<int> int_array1D{};

                std::shared_ptr<anyVar> intvar;
                intvar = std::make_shared<Var<int>>((int)i, varname, current_var_dims,
                                             current_var_atts, vartype, 
                                             int_array1D,
                                             int_array2D);
                vars.push_back(intvar);
            }
            else if (vartype == NC_UINT)
            {
                //////////////////////////////////////////////////////////////////////
                //////////////////////////////////////////////////////////////////////
                uint values[size_y][size_x];
                static size_t start[] = {time_step, 0, 0};
                static size_t count[] = {1, size_y, size_x};
                if ((retval = nc_get_vara_uint(file_id, (int)i, start, count, &values[0][0])))
                    ERR(retval)
                std::vector<std::vector<uint>> uint_array2D;
                for (size_t k_y = 0; k_y < size_y; ++k_y)
                {
                    std::vector<uint> uint_array_row;
                    for (size_t k_x = 0; k_x < size_x; ++k_x)
                    {
                        uint_array_row.push_back(values[k_y][k_x]);
                    }
                    uint_array2D.push_back(uint_array_row);
                }
                //////////////////////////////////////////////////////////////////////
                //////////////////////////////////////////////////////////////////////
                std::vector<uint> uint_array1D{};

                std::shared_ptr<anyVar> uintvar;
                uintvar = std::make_shared<Var<uint>>((int)i, varname, current_var_dims,
                                             current_var_atts, vartype, 
                                             uint_array1D,
                                             uint_array2D);
                vars.push_back(uintvar);
            }
            else
            {
                throw std::invalid_argument("unknown variable type: "+ std::to_string(vartype));
            }
        }
        else if ((nvardims == 1) and time_dep) // ignore time
        {
            //do nothing 
        }
        else if (nvardims == 0)
        {
            std::shared_ptr<anyVar> ndvar;
            std::vector<int> nd_array1D;
            std::vector<std::vector<int>> nd_array2D;
            ndvar = std::make_shared<Var<int>>((int)i, varname, current_var_dims,
                                                current_var_atts, NC_NAT, 
                                                 nd_array1D,
                                                 nd_array2D);
            vars.push_back(ndvar);
        }
        else
        {
            throw std::invalid_argument("can't parce netcdf - unknown variables structure. Currently possible: (x), (y), (t), (y,x), (t,y,x).");
        }
        
    }
}

void Dataset::VarsInfo() const
{
    int i = 0;
    std::cout << "There are " << vars.size() << " variables:" << std::endl;
    for (const std::shared_ptr<anyVar>& item : vars)
    {
        (*item).VarInfo();
    }
}

