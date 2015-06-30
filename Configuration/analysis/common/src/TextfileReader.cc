#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cstdio>
#include <algorithm>

#include "TextfileReader.h"





bool TextfileReader::debug = false;



TextfileReader::TextfileReader():
TextFormatter(),
start_(""),
end_(""),
blindmode_(false),
requirevalues_(true)
{}



void TextfileReader::setStartMarker(std::string d)
{
    d.erase(std::remove_if(d.begin(), d.end(), ::isspace), d.end());
    start_ = d;
}



void TextfileReader::setEndMarker(std::string d)
{
    d.erase(std::remove_if(d.begin(), d.end(), ::isspace), d.end());
    end_ = d;
}



void TextfileReader::readFile(const std::string& filename)
{
    lines_.clear();
    std::ifstream myfile(filename.data(), std::ios::in);
    if(myfile.is_open()){
        tempinfilename_ = filename;
        if(debug){
            std::cout << "TextfileReader::readFile: opened file.\n searching for start marker " << start_
                      << " until end marker " << end_ << "\n"
                      << "delimiter: " << delimiter_ << " comments: " << comment_ << std::endl;

        }
        if(blindmode_){
            std::vector<std::string> record;
            lines_.push_back(record);
        }
        bool started = false;
        bool ended = false;
        while(myfile.good()){
            if(ended) break;
            std::string s;
            
            if(!getline(myfile, s)) break;
            
            if(blindmode_){
                lines_.at(0).push_back(s);
                continue;
            }
            else{
                this->trimcomments(s);
                this->trim(s);
            }
            
            if(s.size() < 1) continue;
            std::string d = s;
            d.erase(std::remove_if(d.begin(), d.end(), ::isspace), d.end());
            // Wait for start marker on line basis
            if(!blindmode_ && !started && start_.size()!=0 && d.find(start_)==std::string::npos) continue;
            else if(!started){
                started = true;
                if(start_.size() > 0){
                    if(debug) std::cout << d << " is start marker, begin reading" << std::endl;
                    // Do not read start marker line itself
                    continue;
                }
            }
            if(started && end_.size()>0 && d.find(end_)!=std::string::npos){
                if(debug) std::cout << d << " is end marker, stop reading" << std::endl;
                ended = true;
                break;
            }
            
            if(debug && !started) std::cout << "found start marker in line" << std::endl;
            
            // Wait for start marker
            //if(end_.size() !=0 && s.compare(end_) == 0) break;

            std::istringstream ss(s);
            std::vector<std::string> record;
            bool noentry = true;
            while(ss){
                std::string s2;
                if(!getline(ss, s2, *delimiter_.data())) break;
                // wait for start marker on word by word basis
                //if(!started && start_.size()!=0 && s2!=start_){
                //    if(debug) std::cout<<s2<<" not yet startmarker"<<std::endl;
                //    continue;
                //}
                //else if(!started){
                //    started = true;
                //    // Do not read start marker line itself
                //    if(start_.size() > 0){
                //        if(debug) std::cout<<s2<<" is start marker, begin reading"<<std::endl;
                //        continue;
                //    }
                //}
                //if(end_.size()>0 && s2==end_){
                //    if(debug) std::cout<<s2<<" is end marker, stop reading"<<std::endl;
                //    ended = true;
                //    break;
                //}
                if(debug) std::cout<<"read \""<<s2<<"\""<<std::endl;
                this->trim(s2);
                if(debug) std::cout<<"trimmed to \""<<s2<<"\""<<std::endl;
                if(s2.size() < 1) continue;
                noentry = false;
                record.push_back(s2);
            }
            if(!ended && !noentry){
                lines_.push_back(record);
                record.clear();
            }
        }
        
        myfile.close();
    }
    else{
        std::cout<<"TextfileReader::readFile: could not read file "<<filename<<std::endl;
        return;
    }
}



std::string TextfileReader::getReJoinedLine(const size_t line)const
{
    if(line >= lines_.size()) throw std::out_of_range("TextfileReader::getReJoinedLine");
    std::string out;
    for(size_t i = 0; i < lines_.at(line).size(); ++i){
        out += lines_.at(line).at(i);
        if(i < lines_.at(line).size()-1) out += delimiter_;
    }
    return out;
}



std::vector<std::string> TextfileReader::getMarkerValues(const std::string& markername)const
{
    std::vector<std::string> out;
    for(size_t i = 0; i < this->nLines(); ++i){
        std::string fullmarker;
        std::vector<std::string> formattedentries;
        if(delimiter_==" " || delimiter_=="-" || delimiter_=="]" ||delimiter_=="["){
            fullmarker = this->getReJoinedLine(i);
            TextFormatter tf;
            fullmarker = this->getReJoinedLine(i);
            tf.setComment(comment_);
            tf.setDelimiter(",");
            tf.setTrim(trim_);
            formattedentries = tf.getFormatted(fullmarker);
        }
        else{
            formattedentries = lines_.at(i);
        }
        for(size_t entr = 0; entr < formattedentries.size(); ++entr){
            TextFormatter tfentr;
            tfentr.setDelimiter("-");
            tfentr.setComment(comment_);
            tfentr.setTrim(trim_);
            std::vector<std::string> varsandvals = tfentr.getFormatted(formattedentries.at(entr));
            
            if(varsandvals.size() < 1) continue;
            // no marker
            if(varsandvals.at(0).find("[") ==  std::string::npos) continue;
            tfentr.setTrim("[");
            tfentr.ltrim(varsandvals.at(0));
            tfentr.setTrim(" ");
            tfentr.trim(varsandvals.at(0));
            if(varsandvals.at(0) == markername){
                tfentr.setTrim("[]");
                tfentr.trim(formattedentries.at(entr));
                tfentr.setTrim(markername);
                tfentr.trim(formattedentries.at(entr));
                tfentr.setTrim("- ");
                tfentr.trim(formattedentries.at(entr));
                out.push_back(formattedentries.at(entr));
            }
        }
    }
    return out;
}



std::vector<std::string> TextfileReader::getData(const size_t& line)const
{
    if(line >= lines_.size()){
        std::cout<<"TextfileReader::getData: line out of range"<<std::endl;
        return std::vector<std::string>();
    }
    return lines_.at(line);
}



std::string TextfileReader::getValueString(const std::string& val, const bool checkdoubles)
{
    std::string out = "";
    size_t count = 0;
    // search all entries
    for(size_t line = 0; line < this->nLines(); ++line){
        std::string thisline;
        TextFormatter tf;
        std::vector<std::string> formattedentries;
        if(delimiter_ != ","){
            thisline = this->getReJoinedLine(line);
            tf.setComment(comment_);
            tf.setDelimiter(",");
            tf.setTrim(trim_);
            formattedentries = tf.getFormatted(thisline);
        }
        else{
            formattedentries = lines_.at(line);
        }
        
        //split according to format, allow for commas as separators and spaces in names
        for(size_t entr = 0; entr < formattedentries.size(); ++entr){
            TextFormatter tfentr;
            tfentr.setDelimiter("=");
            tfentr.setComment(comment_);
            tfentr.setTrim(trim_);
            std::vector<std::string> varsandvals = tfentr.getFormatted(formattedentries.at(entr));
            if(varsandvals.size() < 2) continue;
            
            if(varsandvals.at(0) == val){ //found
                if(checkdoubles && count>0){//except
                    std::cout<<"TextfileReader::getValue: value "<<val<<" defined twice! Source of errors!"<<std::endl;
                    throw std::logic_error("TextfileReader::getValue: value defined twice! Source of errors!");
                }
                out = varsandvals.at(1);
                ++count;
            }
        }
    }
    return out;
}



std::string TextfileReader::dumpFormattedToTmp()const
{
    char buffer [256] = "/tmp/tmpfileXXXXXX";
    mkstemp(buffer);
    std::string filename(buffer);
    std::ofstream ofs(filename, std::ofstream::out);
    for(size_t i = 0; i < this->nLines(); ++i){
        for(size_t entr = 0; entr < lines_.at(i).size(); ++entr){
            ofs << lines_.at(i).at(entr) << delimiter_;
        }
        ofs << std::endl;
    }
    ofs.close();
    return filename;
}




