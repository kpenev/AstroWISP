#include "FitsHeader.h"

namespace IO {

    static void trim(std::string &s,
                     const std::string &trim_characters=" \f\n\r\t\v")
    {
        s.erase(0, s.find_first_not_of(trim_characters));
        s.erase(s.find_last_not_of(trim_characters)+1);
    }

    void FitsHeader::read(fitsfile *fptr)
    {
        int num_records,
            fits_status = 0;
        fits_get_hdrspace(fptr, &num_records, NULL, &fits_status);
        if(fits_status) 
            throw Error::Fits(
                "Failed to count the number of keyword in fits header."
            );
        std::string prev_keyword;
        for(int record = 1; record <= num_records; ++record) {
            char keyword[80],
                 value[80],
                 comment[80];
            fits_read_keyn(fptr,
                           record,
                           keyword,
                           value,
                           comment,
                           &fits_status);
            if(fits_status) {
                std::ostringstream msg;
                msg << "Failed to read keyword number " << record
                    << " from fits header.";
                throw Error::Fits(msg.str().c_str());
            }
            std::string keyword_str(keyword);
            std::string clean_val;
            if(keyword_str == "COMMENT" || keyword_str == "CONTINUE") {
                clean_val = comment;
                trim(clean_val);
                if(clean_val[0] == '=') clean_val.erase(clean_val.begin());
            } else clean_val = value;
            trim(clean_val);
            if(keyword_str == "CONTINUE") keyword_str = prev_keyword;
            else prev_keyword = keyword_str;
            if(clean_val[0] == '\'') clean_val.erase(clean_val.begin());
            if(clean_val[0] == '>') clean_val.erase(clean_val.begin());
            if(clean_val[clean_val.size()-1] == '\'')
                clean_val.erase(clean_val.size() - 1);
            trim(clean_val);
            std::string *dest = &__values[keyword_str];
            if(dest->size() == 0) {
                (*dest) = clean_val;
                __keywords.push_back(keyword_str);
                __comments[keyword_str] = comment;
            } else if((*dest)[dest->size() - 1] == '\\') {
                dest->erase(dest->size() - 1);
                dest->append(clean_val);
            } else {
                dest->push_back(keyword_str == "MASKINFO" ? ' ' : '\n');
                dest->append(clean_val);
            }
        }
    }

} //End IO namespace.
