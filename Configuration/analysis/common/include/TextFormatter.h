#ifndef TextFormatter_h
#define TextFormatter_h

#include <string>
#include <vector>
#include <iostream>
#include <sstream>





/// Class for formatting text with various functionalities
class TextFormatter{
    
public:
    
    /// Constructor
    TextFormatter();
    
    /**
     * sets all chars that are used for trimming
     */
    void setTrim(const std::string& tr){trim_ = tr;}
    
    /**
     * defines one char that will serve as comment indicator
     */
    void setComment(const std::string& c){comment_ = c;}
    
    /**
     * sets a delimiter for individual entries per line
     */
    void setDelimiter(const std::string& d){delimiter_ = d;}
    
    /**
     * cuts on left side of the string
     * example:
     *     setTrim(#)
     *     string="#######this is a string####"
     *     ltrim(string)
     *     now: string="this is a string####"
     */
    std::string& ltrim(std::string& str)const;
    
    /**
     * cuts on right side of the string
     * example:
     *     setTrim(#)
     *     string="#######this is a string####"
     *     rtrim(string)
     *     now: string="#######this is a string"
     */
    std::string& rtrim(std::string& str)const;
    
    /// trim from both ends
    std::string& trim(std::string& str)const;
    
    /**
     * cuts everything that follows the comment marker (one char)
     * example:
     *     setComment("&")
     *     string="this is a string &with a comment"
     *     trimcomments(string)
     *     now: string=="this is a string "
     */
    std::string& trimcomments(std::string& str)const;
    
    
    /**
     * Template for returning formatted text fragments according to the format rules defined previously
     */
    template<class T>
    std::vector<T> getFormatted(const std::string& in)const;
    
    /**
     * Return formatted text fragments as std::string according to the format rules defined previously
     */
    std::vector<std::string> getFormatted(const std::string& in)const{return this->getFormatted<std::string>(in);}
    
    /**
     * returns file name with extension but without path
     */
    static std::string getFilename(const std::string& pathtofile);
    
    /**
     * returns file extension or empty string if none
     */
    static std::string getFileExtension(const std::string& pathtofile);
    
    /**
     * strips file of extension
     */
    static std::string stripFileExtension(const std::string& pathtofile);
    
    /**
     * strips file of dir
     */
    static std::string stripFileDir(const std::string& pathtofile);
    
    /**
     * gets file directory
     */
    static std::string getFileDir(const std::string& pathtofile);
    
    /**
     * adds a suffix to file name but keeps extension, e.g. path/bla.txt is transformed to
     * path/bla_suffix.txt
     */
    static std::string addFilenameSuffix(const std::string& pathtofile, const std::string& suffix);
    
    /// Replace all special characters by underscores
    static std::string makeCompatibleFileName(const std::string& in);
    
    /**
     * Switch for more output
     */
    static bool debug;
    
    
    
protected:
    
    /// Characters used for trimming, can be more than one, e.g. "\t\n "
    std::string trim_;
    
    /// Char that will serve as comment indicator
    /// All text following that char will be ignored until the next line starts
    std::string comment_;
    
    /// Delimiter for individual entries per line (e.g. "," or " ")
    std::string delimiter_;
};





template<class T>
inline std::vector<T> TextFormatter::getFormatted(const std::string& in)const
{
    std::vector<T> out;
    std::string s = in;
    this->trimcomments(s);
    this->trim(s);
    std::istringstream ss(s);
    while(ss){
        std::string s2;
        if(!getline(ss, s2, *delimiter_.data())) break;
        if(debug) std::cout<<"got \""<<s2<<"\""<<std::endl;
        this->trim(s2);
        if(debug) std::cout<<"trimmed to \""<<s2<<"\""<<std::endl;
        if(s2.length() > 0) out.push_back((T)s2);
    }
    return out;
}



#endif




