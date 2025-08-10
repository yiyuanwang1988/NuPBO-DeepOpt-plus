#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>

using namespace std;

void translator(string opb_file_name, string output_file_name)
{
    const int max_var = 50000000;   //约束结构体可以存储的最大变量数
    struct constraint_type        //约束结构类型
    {
        int coeff[max_var];
        int var[max_var];
        string symbol;
        int degree;
    } constraint, equal_constraint;

    long long num_vars, num_clauses, num_softclauses, num_hardclauses, num_constraints;
    long long top_clause_weight;

    ifstream opb_ifile(opb_file_name.c_str());      //ifstream是从硬盘到内存，c_str返回当前字符串的首字符地址
    ofstream wcnf_ofile(output_file_name.c_str());  //ofstream是从内存到硬盘

    stringstream ss;             //strstream类：字符串流的输入输出操作
    stringstream ss2;

    stringstream mss;
    mss.clear();

    string line, firstline, other;

    int i = 0;              //i是循环变量
    int sum = 0;            //sum是系数之和
    int negsum = 0;         //negsum是负系数之和
    int equal_cnt = 0;      //equal_cnt是等号情况的数目

    getline(opb_ifile, line);   //p cnf 1338 2249，从输入流接收一行字符串，并放入字符串line中
    ss.clear();                 //多次使用stringstream，需要清空
    ss.str(line);               //把line中的字符串存入字符串流中 
    ss.seekg(0, ios::beg);      //seekg(offset, place)设置输入文件流的文件流指针位置，ios::beg文件的开始位置
    ss >> other;                //过滤掉 *
    ss >> other;                //过滤掉 #variable=
    ss >> num_vars;             //获取变量数
    ss >> other;                //过滤掉 #constraint=
    ss >> num_hardclauses;      //获取硬子句数目（即opb文件中的约束数目）

    //处理权重、子句数目和软子句
    //第1行：min: +1 x1 +2 x2 +4 x3 +8 x4 +6 x38 ;

    string coeffl, varl;
    getline(opb_ifile, line);
    while(line[0]=='*') {getline(opb_ifile, line);} //过滤掉注释，*开头

    //计算权重和软子句个数;
    ss.clear();
    ss.str(line);           
    ss.seekg(0, ios::beg); 
    ss >> coeffl;                //过滤掉第1行中的“min:“
    ss >> coeffl >> varl;         //coeff：+1，var：x1
    long long icoeff, ivar;
    num_softclauses = 0;        //软子句数目为0
    top_clause_weight = 0;      //最高子句权重为0
    while (coeffl != ";")        //未到一个子句末尾，；为子句结束符
    {
        ss2.clear();
        ss2.str(coeffl);         //+1
        ss2.seekg(0, ios::beg);
        ss2 >> icoeff;          //+1
        ss2.clear();
        ss2.str(varl.substr(1)); //过滤掉“x”；substr(start, length)：返回一个从指定位置开始，并具有指定长度的子字符串，参数start必选，字符串中第一个字符的索引为0
        ss2.seekg(0, ios::beg);
        ss2 >> ivar;            //1
        num_softclauses++;      //软子句数目加1
        if(icoeff>0)                          //+1 > 0
        	top_clause_weight += icoeff;      //top_clause_weight = 0 + (+1)
        else 
            top_clause_weight += (-icoeff);   //top_clause_weight = 0 - (-1)
        ss >> coeffl >> varl;     //coeff：+2，var：x2
    }

    firstline = line;           //保留软子句对应的第1行
    
    //处理硬子句
    //第2行
    while(getline(opb_ifile, line))
    {
        /* 扫面一行，存入constraint */
        ss.clear();
        ss.str(line);
        ss.seekg(0, ios::beg);
        ss >> coeffl >> varl;           //coeff：+1，var：x1
        
        i = 0;
        while (coeffl != ">=" && coeffl != ">" && coeffl != "<=" && coeffl != "<" && coeffl != "=") 
        {            
            ss2.clear();
            ss2.str(coeffl);
            ss2.seekg(0, ios::beg);
            ss2 >> icoeff;            //将字符串coeff转化为整型icoeff
            
            constraint.coeff[i] = icoeff;

            ss2.clear();
            ss2.str(varl.substr(1));
            ss2.seekg(0, ios::beg);
            ss2 >> ivar;              //将字符串var转化为整型ivar

            constraint.var[i] = ivar;
            
            i++;

            ss >> coeffl >> varl;
        }
        
        constraint.coeff[i] = 0;      //系数数组末尾标志0
        constraint.var[i] = 0;        //变量数组末尾标志0

        constraint.symbol = coeffl;    //最后一个coeff存放的是symbol

        ss2.clear();
        ss2.str(varl);
        ss2.seekg(0, ios::beg);
        ss2 >> ivar;
        constraint.degree = ivar;     //最后一个var存放的是degree

        /* 正则化第一步：符号变换 */

        if (constraint.symbol == ">")    //处理">"的情况
        {
            constraint.symbol = ">=";
            constraint.degree += 1;
        } 
        else if (constraint.symbol == "<=")    //处理"<="的情况
        {
            i = 0;
            sum = 0;
            while (constraint.coeff[i] != 0)
            {
                sum += constraint.coeff[i];
                i++;
            } 
            
            i = 0;
            while (constraint.coeff[i] != 0)
            {
                constraint.var[i] = -constraint.var[i];
                i++;
            }

            constraint.symbol = ">=";

            constraint.degree = sum - constraint.degree;
        }
        else if (constraint.symbol == "<")      //处理"<"的情况
        {
            i = 0;
            sum = 0;
            while (constraint.coeff[i] != 0)
            {
                sum += constraint.coeff[i];
                i++;
            }
            
            i = 0;
            while (constraint.coeff[i] != 0)
            {
                constraint.var[i] = -constraint.var[i];
                i++;
            }

            constraint.symbol = ">=";

            constraint.degree = sum - constraint.degree + 1;
        }
        else if (constraint.symbol == "=")       //处理"="的情况
        {
            equal_cnt++;
            //转换的第1个约束 无变化，即constraint

            //转换的第2个约束 构造equal_constraint
            i = 0;
            sum = 0;
            while (constraint.coeff[i] != 0)
            {
                equal_constraint.coeff[i] = constraint.coeff[i];
                sum += constraint.coeff[i];
                i++;
            }
            equal_constraint.coeff[i] = 0;
          
            i = 0;
            while (constraint.coeff[i] != 0)
            {
                equal_constraint.var[i] = -constraint.var[i];
                i++;
            }
            equal_constraint.var[i] = 0;

            equal_constraint.symbol = ">=";

            equal_constraint.degree = sum - constraint.degree; 
        }

        /* 正则化第二步：消除负系数 */
        negsum = 0;    //消除constraint负系数
        i = 0;
        while (constraint.coeff[i] != 0)
        {
            if (constraint.coeff[i] < 0)
            {
                negsum += -constraint.coeff[i];
                constraint.coeff[i] = -constraint.coeff[i];
                constraint.var[i] = -constraint.var[i];
            }   
            i++;
        }
        
        constraint.degree += negsum;

        if (constraint.symbol == "=")
        {
            negsum = 0;    //消除equal_constraint负系数
            i = 0;
            while (equal_constraint.coeff[i] != 0)
            {
                if (equal_constraint.coeff[i] < 0)
                {
                    negsum += -equal_constraint.coeff[i];
                    equal_constraint.coeff[i] = -equal_constraint.coeff[i];
                    equal_constraint.var[i] = -equal_constraint.var[i];
                }   
                i++;
            }

            equal_constraint.degree += negsum; 
        }
          
        /* 写入约束 */   
        mss << top_clause_weight+1 << " ";    //写入权重；

        mss << constraint.degree << " ";    //写入度； 

        i = 0;
        while (constraint.coeff[i] != 0)           //写入系数和变量
        {
            mss << constraint.coeff[i] << " " << constraint.var[i] << " ";
            i++;
        }

        mss << "0" << endl;                 //写入结束符

        if (constraint.symbol == "=")              //处理“=”情况，写入转换的第2个约束
        {
            mss << top_clause_weight+1 << " ";

            mss << equal_constraint.degree << " ";  

            i = 0;
            while (equal_constraint.coeff[i] != 0) 
            {
                mss << equal_constraint.coeff[i] << " " << equal_constraint.var[i] << " ";
                i++;
            }

            mss << "0" << endl;
        }                
    }

    //处理软子句
    //第1行
    ss.clear();
    ss.str(firstline);
    ss.seekg(0, ios::beg);
    ss >> coeffl;                //过滤掉第1行中的“min:“
    ss >> coeffl >> varl;         //coeff：+1，var：x1
    while (coeffl != ";")
    {
        ss2.clear();
        ss2.str(coeffl);         //+1
        ss2.seekg(0, ios::beg);
        ss2 >> icoeff;          //+1
        ss2.clear();
        ss2.str(varl.substr(1)); //1
        ss2.seekg(0, ios::beg);
        ss2 >> ivar;            //1
        if (icoeff < 0)
        {
            mss << -icoeff << " 1" << " 1" << " " << ivar << " 0" << endl; //-1 1 0
        }
        else
        {
            mss << icoeff << " 1" << " 1" << " " << -ivar << " 0" << endl; //1 -1 0
        }
        ss >> coeffl >> varl;     //coeff：+2，var：x2
    }

    wcnf_ofile << "p wcnf " << num_vars << " " << num_hardclauses + num_softclauses + equal_cnt << " " << top_clause_weight+1 << endl;  //p wcnf 1138 2249+1 1+1
    wcnf_ofile << mss.str();

    opb_ifile.close();
    wcnf_ofile.close();
}

int main(int argc, char *argv[])
{
    translator(argv[1], argv[2]);
}
