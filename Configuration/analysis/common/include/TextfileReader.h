#ifndef TextTextfileReader_h
#define TextTextfileReader_h

#include <string>
#include <vector>
#include <sstream>
#include <stdexcept>
#include <iostream>

#include <TString.h>

#include "TextFormatter.h"






/// Class for reading in formatted text files
class TextfileReader : public TextFormatter{
    
public:
    
    /// Constructor
    TextfileReader();
    
    /// Destructor
    ~TextfileReader(){this->clear();}
    
    /**
     * Set start marker to be searched in file
     * Ignores any white spaces
     */
    void setStartMarker(std::string d);
    
    /**
     * Set end marker to be searched in file
     * Ignores any white spaces
     */
    void setEndMarker(std::string d);
    
    /// Switch to blind reading mode, ignoring all delimiters, comments and trimmers
    void setBlindMode(bool blind){blindmode_ = blind;}
    
    /// Set if for parameters read in values are required, or if empty argument is allowed
    void setRequireValues(bool req){requirevalues_ = req;}
    
    /// Read in text file with given name
    void readFile(const std::string& filename);
    
    /// Check if vector of all lines is empty
    bool isEmpty()const{return lines_.size() < 1;}
    
    /// Get number of lines
    size_t nLines()const{return lines_.size();}
    
    /// Get number of entries for specified line
    size_t nEntries(const size_t& line)const{return this->getData(line).size();}
    
    /// Get all arguments in specified line
    std::vector<std::string> getData(const size_t& line)const;
    
    /// Get value of a certain entry in a certain line
    template<class T>
    T getData(const size_t& line, const size_t& entry)const;
    
    /// Join all entries from one line into one string, separated by delimiter, and return it
    std::string getReJoinedLine(const size_t line)const;
    
    /**
     * returns additional information connected to a marker formatted:
     *
     * [markername - markervalue] // additional spaces are allowed everywhere
     *
     */
    std::vector<std::string> getMarkerValues(const std::string& markername)const;
    
    /**
     * if file has entries like somevariable=blabla
     * getValue("somevariable") will return "blabla"
     * if value is definied several times, an exception is thrown
     * less performance but safer
     * can be switched off by bool
     * then last value is returned
     *
     * if value is not found empty string will be returned
     *
     *
     * fixed format: commas as separators (if any)
     */
    std::string getValueString(const std::string& str, const bool checkdoubles =true);
    
    /// Get value of certain parameter
    template<class T>
    T getValue(const std::string& str);
    
    /// Get value of certain parameter, or if requirevalues_==false use specified default value
    template<class T>
    T getValue(const std::string& str, T def_val);
    
    /// Clear the texts which were read in
    void clear(){lines_.clear();}
    
    /**
     * Switch for more output
     */
    static bool debug;
    
    /**
     * Write formatted text to file and store in /tmp
     * Returns path to file
     */
    std::string dumpFormattedToTmp()const;
    
    
    
private:
    
    // make inherited private
    std::vector<std::string> getFormatted(const std::string&)const{return std::vector<std::string>();}
    
    
    /// Start marker for identifying part in text file to be read
    std::string start_;
    
    /// End marker for identifying part in text file to be read
    std::string end_;
    
    /// Vector holding for each line a vector with each value in this line
    std::vector<std::vector<std::string> > lines_;
    
    /// Read file in blind mode, i.e. read everything ignoring meaning of all delimiters, comments, trimmers
    /// Everything is written to the first line in lines_
    bool blindmode_;
    
    /// Whether values are required when reading a parameter, or whether defaults are allowed
    bool requirevalues_;
    
    /// Name of file to be read
    std::string tempinfilename_;
};





template<class T>
inline T TextfileReader::getData(const size_t& line, const size_t& entry)const
{
    if(line >= lines_.size()){
        throw std::out_of_range("TextfileReader::getData: line out of range");
    }
    if(entry >= lines_.at(line).size()){
        throw std::out_of_range("TextfileReader::getData: entry out of range");
    }
    T out;
    std::stringstream ss(this->getData(line).at(entry));
    ss >> out;
    return out;
}



template<class T>
inline T TextfileReader::getValue(const std::string& str)
{
    T out;
    std::string s(this->getValueString(str, true));
    if(s.length() < 1){
        std::string exstring = "TextfileReader::getValue: value for " + str + " not found in " + tempinfilename_;
        throw std::runtime_error(exstring);
    }
    std::stringstream ss(s);
    ss >> out;
    return out;
}



template<class T>
inline T TextfileReader::getValue(const std::string& str, T def_val)
{
    if(requirevalues_) return this->getValue<T>(str);
    T out;
    std::string s(this->getValueString(str, true));
    if(s.length() > 0){
        std::stringstream ss(s);
        ss >> out;
        return out;
    }
    else{
        return def_val;
    }
}



template<>
inline bool TextfileReader::getData<bool>(const size_t& line, const size_t& entry)const
{
    bool out;
    std::stringstream ss(this->getData(line).at(entry));
    ss >> std::boolalpha >> out;
    return out;
}



template<>
inline bool TextfileReader::getValue<bool>(const std::string& str)
{
    bool out;
    std::string s(this->getValueString(str, true));
    if(s.length() < 1){
        std::string exstring = "TextfileReader::getValue: value for " + str + " not found in " + tempinfilename_;
        throw std::runtime_error(exstring);
    }
    std::stringstream ss(s);
    ss >> std::boolalpha >> out;
    return out;
}



template<>
inline bool TextfileReader::getValue<bool>(const std::string& str, bool def_val)
{
    if(requirevalues_) return this->getValue<bool>(str);
    bool out;
    std::string s(this->getValueString(str, true));
    if(s.length() > 0){
        std::stringstream ss(s);
        ss >> std::boolalpha >> out;
        return out;
    }
    else{
        return def_val;
    }
}



template<>
inline std::string TextfileReader::getData<std::string>(const size_t& line,const size_t& entry)const
{
    return (this->getData(line).at(entry));
}



template<>
inline std::string TextfileReader::getValue<std::string>(const std::string& str)
{
    std::string s(this->getValueString(str));
    if(s.length() < 1){
        std::string exstring = "TextfileReader::getValue: value for " + str + " not found in " + tempinfilename_;
        throw std::runtime_error(exstring);
    }
    return s;
}



template<>
inline std::string TextfileReader::getValue<std::string>(const std::string& str, std::string def_val)
{
    if(requirevalues_) return this->getValue<std::string>(str);
    std::string s(this->getValueString(str));
    if(s.length() < 1) return def_val;
    return s;
}



template<>
inline TString TextfileReader::getData<TString>(const size_t& line,const size_t& entry)const
{
    return (TString)(this->getData(line).at(entry));
}



template<>
inline TString TextfileReader::getValue<TString>(const std::string& str)
{
    std::string s(this->getValueString(str));
    if(s.length() < 1){
        std::string exstring = "TextfileReader::getValue: value for " + str + " not found in " + tempinfilename_;
        throw std::runtime_error(exstring);
    }
    return (TString)s;
}



template<>
inline TString TextfileReader::getValue<TString>(const std::string& str, TString def_val)
{
    if(requirevalues_) return this->getValue<TString>(str);
    std::string s(this->getValueString(str));
    if(s.length() < 1) return def_val;
    return (TString)s;
}





#endif



