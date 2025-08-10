#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>

using namespace std;

void convertor(string opb_file_name, string cnf_file_name, string output_file_name)
{
    long long num_vars, num_clauses, num_softclauses, num_hardclauses, num_constraints;
    long long top_clause_weight;

    ifstream opb_ifile(opb_file_name.c_str());      //ifstream是从硬盘到内存
    ifstream cnf_ifile(cnf_file_name.c_str());      //c_str返回当前字符串的首字符地址
    ofstream wcnf_ofile(output_file_name.c_str());  //ofstream是从内存到硬盘

    stringstream ss; //strstream类：字符串流的输入输出操作
    string line, line2, line3;

    getline(cnf_ifile, line); //p cnf 1338 2249，从输入流接收一行字符串，并放入字符串line中
    ss.clear();   //多次使用stringstream，需要清空
    ss.str(line); //把line中的字符串存入字符串流中 
    ss.seekg(0, ios::beg); //seekg(offset, place)设置输入文件流的文件流指针位置，ios::beg文件的开始位置
    ss >> line2 >> line3 >> num_vars >> num_hardclauses; //line2：p，line3：cnf， num_vars：1338，num_hardclauses：2249

    //处理软子句

    string coeff, var;
    getline(opb_ifile, line);
    while(line[0]=='*') {getline(opb_ifile, line);} //过滤注释，*开头
    //getline(opb_ifile, line);
    //第1行：min: +1 x1 +2 x2 +4 x3 +8 x4 +6 x38 ;
    ss.clear();
    ss.str(line);           
    ss.seekg(0, ios::beg); 
    ss >> coeff;           //过滤第1行中的“min:“，min: +1 x1 +2 x2 +4 x3 +8 x4 +6 x38 ;
    ss >> coeff >> var;    //coeff：+1，var：x1
    long long icoeff, ivar;
    num_softclauses = 0;   //软子句数为0
    top_clause_weight = 0; //最高子句权重为0
    while (coeff != ";")   //未到一个子句末尾，；为子句结束符
    {
        stringstream ss2;
        ss2.clear();
        ss2.str(coeff);    //+1
        ss2.seekg(0, ios::beg);
        ss2 >> icoeff;     //+1
        ss2.clear();
        ss2.str(var.substr(1)); //过滤掉“x”，substr(start, length)：返回一个从指定位置开始，并具有指定长度的子字符串，参数start必选，字符串中第一个字符的索引为0
        ss2.seekg(0, ios::beg);
        ss2 >> ivar;       //1
        num_softclauses++; //软子句个数加1
        if(icoeff>0)                          //+1 > 0
        	top_clause_weight += icoeff;      //top_clause_weight = 0 + (+1)
        else top_clause_weight += (-icoeff);  //top_clause_weight = 0 - (-1)
        ss >> coeff >> var;     //coeff：+2，var：x2
    }
    wcnf_ofile << "p wcnf " << num_vars << " " << num_hardclauses + num_softclauses << " " << ++top_clause_weight << endl;  //p wcnf 1138 2249+1 1+1

    //第2+行：+1 x1 +2 x2 +4 x3 +8 x4 +1 x5 +2 x6 +4 x7 +8 x8 -2 x14 -4 x15 -8 x16 -16 x17 -16 x37 >= -16 ; ...
    ss.clear();
    ss.str(line);
    ss.seekg(0, ios::beg);
    ss >> coeff;           //???
    ss >> coeff >> var;    //coeff：+1，var：x1
    while (coeff != ";")
    {
        stringstream ss2;
        ss2.clear();
        ss2.str(coeff);    //+1
        ss2.seekg(0, ios::beg);
        ss2 >> icoeff;     //+1
        ss2.clear();
        ss2.str(var.substr(1)); //1
        ss2.seekg(0, ios::beg);
        ss2 >> ivar;       //1
        if (icoeff < 0)
        {
            wcnf_ofile << -icoeff << " " << ivar << " 0" << endl; //-1 1 0
        }
        else
        {
            wcnf_ofile << icoeff << " " << -ivar << " 0" << endl; //1 -1 0
        }
        ss >> coeff >> var; //coeff：+2，var：x2
    }
    opb_ifile.close();

    //处理硬子句

    for (int i = 0; i < num_hardclauses; ++i)
    {
        getline(cnf_ifile, line);
        wcnf_ofile << top_clause_weight << " " << line << endl; //1+2+4+8+6+1=22 39 0
    }
    cnf_ifile.close();
    wcnf_ofile.close();
}

int main(int argc, char *argv[])
{
    convertor(argv[1], argv[2], argv[3]);
}
