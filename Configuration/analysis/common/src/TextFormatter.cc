#include <iostream>
#include <sstream>
#include <stdexcept>
#include <algorithm>

#include "TextFormatter.h"





bool TextFormatter::debug = false;



TextFormatter::TextFormatter():
trim_(" \t"),
comment_("#"),
delimiter_(",")
{}



std::string& TextFormatter::rtrim(std::string& str)const
{
    if(trim_.length() < 1) return str;
    size_t endpos = str.find_last_not_of(trim_);
    if(std::string::npos != endpos) str = str.substr(0, endpos+1);
    else str.clear(); //only contains trim chars
    return str;
}



std::string& TextFormatter::ltrim(std::string& str)const
{
    if(trim_.length() < 1) return str;
    size_t startpos = str.find_first_not_of(trim_);
    if(std::string::npos != startpos) str = str.substr(startpos );
    else str.clear(); //only contains trim chars
    return str;
}



std::string& TextFormatter::trimcomments(std::string& str)const
{
    if(comment_.length() < 1) return str;
    if(str.length()<2 && str==comment_) str = "";
    size_t endpos = str.find(comment_);
    if(std::string::npos != endpos) str = str.substr( 0, endpos);
    return str;
}



std::string& TextFormatter::trim(std::string& s)const
{
    this->ltrim(this->rtrim(s));
    return s;
}





// ------------------------- static member functions -------------------------


std::string TextFormatter::getFilename(const std::string& pathtofile)
{
    std::string out;
    std::string s = pathtofile;
    std::istringstream ss(s);
    while(ss){
        std::string s2;
        if(!getline(ss, s2, *"/")) break;
        if(debug) std::cout<<"got \""<<s2<<"\""<<std::endl;
        if(s2.size() > 0) out = s2;
    } // just keep last
    return out;
}



std::string TextFormatter::getFileExtension(const std::string& pathtofile)
{
    std::string out = getFilename(pathtofile);
    std::string s = out;
    std::istringstream ss(s);
    while(ss){
        std::string s2;
        if(!getline( ss, s2, *".")) break;
        if(debug) std::cout<<"got \""<<s2<<"\""<<std::endl;
        if(s2.size() > 0) out = s2;
    } // just keep last
    //in case no extension
    if(out == s) return "";
    return out;
}



std::string TextFormatter::stripFileExtension(const std::string& pathtofile)
{
    std::string out = "";
    std::string s = getFilename(pathtofile);
    std::istringstream ss(s);
    size_t endpos = pathtofile.find_last_of(".");
    if(std::string::npos != endpos) out = pathtofile.substr(0, endpos);
    return out;
}



std::string TextFormatter::stripFileDir(const std::string& pathtofile)
{
    std::string out = "";
    std::string s = getFilename(pathtofile);
    std::istringstream ss(s);
    size_t endpos = pathtofile.find_last_of("/");
    if(std::string::npos != endpos) out = pathtofile.substr(endpos+1, pathtofile.length());
    return out;
}



std::string TextFormatter::getFileDir(const std::string& pathtofile){
    std::string str = pathtofile;
    size_t endpos = str.find_last_of("/");
    if(std::string::npos != endpos) str = str.substr(0, endpos+1);
    return str;
}



std::string TextFormatter::addFilenameSuffix(const std::string& pathtofile, const std::string& suffix)
{
    std::string onlyname = stripFileExtension(pathtofile);
    std::string extension = getFileExtension(pathtofile);
    return onlyname + suffix + "." + extension;
}



std::string TextFormatter::makeCompatibleFileName(const std::string& in)
{
    std::string out = in;
    std::replace(out.begin(), out.end(), '#', '_');
    std::replace(out.begin(), out.end(), '/', '_');
    std::replace(out.begin(), out.end(), '{', '_');
    std::replace(out.begin(), out.end(), '}', '_');
    std::replace(out.begin(), out.end(), ' ', '_');
    std::replace(out.begin(), out.end(), '\\', '_');
    return out;
}



