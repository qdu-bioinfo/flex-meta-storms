// Updated at Sept 30, 2022
// Updated by Mingqian Zhang
// Bioinformatics Group, College of Computer Science & Technology, Qingdao University
// Code by: Mingqian Zhang, Xiaoquan Su
// Flex Meta-Storms version 1.0

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <vector>

#include "fms_func.h"
#include "fms_db.h"
#include "hash.h"
#include "fms_version.h"
#include "table_format.h"
#include "otu_parser.h"
#include "dist.h"
#include "comp.h"

#ifndef FMS_COMP_H
#define FMS_COMP_H

#define REG_SIZE 70

#define MIN_DIST 0.00001

using namespace std;

class _FMS_Comp_Tree : public _Comp_Tree{
      
      public:
             _FMS_Comp_Tree(){                           
                           LeafN = 0;
                           OrderN = 0;                 
                           Init();
                          }
    
            _FMS_Comp_Tree(char db){
                        Database.Set_DB(db);  
						
                        LeafN = 0;
                        OrderN = 0;
						Init();
                        }
    
             int Load_FMS_abd(_Table_Format * table, float * Abd, int sample, bool is_cp_correct); //Load by table_format
             float Calc_FMS_sim(float * Abd_1, float * Abd_2);
             int Get_LeafN(){
                 return LeafN;
                 }
             
      private:
              _PMDB Database;
    
              vector <float> Dist_1;
              vector <float> Dist_2;
              
              vector <int> Order_1;
              vector <int> Order_2;
              vector <int> Order_d;
              
              vector <string> Id; 
              
              hash_map <string, float, std_string_hash> Cp_number;
              hash_map <string, int, std_string_hash> Id_hash;
    
              void Init();
              int Load_id();
              int Load_order();   
              
              int LeafN;
              int OrderN;                         
              };

void _FMS_Comp_Tree::Init(){    
                    
     LeafN = Load_id();
     //load id hash
     for (int i = 0; i < LeafN; i ++)
        Id_hash[Id[i]] = i;
        
     OrderN = 0;
     //load tree  
     if (Database.Get_Is_Tree())      
        OrderN = Load_order();
    
     //load cp number         
     if (Database.Get_Is_Cp())
        Database.Load_Copy_Number(Cp_number);
     }

int _FMS_Comp_Tree::Load_id(){
	
     ifstream infile(Database.Get_Tree_Id().c_str(), ifstream::in);
     if (!infile){
                 cerr << "Error: Cannot open file : " << Database.Get_Tree_Id() << endl;
                 return 0;
                 }
     
     string buffer;
     int count = 0;
     while(getline(infile, buffer)){
                           if (buffer.size() == 0) continue;
                           Id.push_back(buffer);
                           count ++;
                           }
     
     infile.close();
     infile.clear();
     return count;
     }


int _FMS_Comp_Tree::Load_order(){
	
    ifstream infile(Database.Get_Tree_Order().c_str(), ifstream::in);
    if (!infile){
                 cerr << "Error: Cannot open file : " << Database.Get_Tree_Order() << endl;
                 return 0;
                 }
    
     string buffer;
     int count = 0;
     while(getline(infile, buffer)){                           
                           if(buffer.size() == 0) continue;
                           stringstream strin(buffer);
                           int order_1 = 0;
                           int order_2 = 0;
                           int order_d = 0;
                           float dist_1 = 0;
                           float dist_2 = 0;
                           strin >> order_1 >> dist_1 >> order_2 >> dist_2 >> order_d;
                           Order_1.push_back(order_1);
                           Order_2.push_back(order_2);
                           Order_d.push_back(order_d);
                           Dist_1.push_back(dist_1);
                           Dist_2.push_back(dist_2);
                           count ++;                           
                           }
       
    infile.close();
    infile.clear(); 
    
    return count;
    }

int _FMS_Comp_Tree::Load_FMS_abd(_Table_Format * table, float *Abd, int sample, bool is_cp_correct){
    
    memset(Abd, 0, LeafN * sizeof(float));

    vector <string> otus = table->Get_Feature_Names();
    vector <float> abds = table->Get_Abd(sample);
    
    hash_map<string, float, std_string_hash> otu_abd;
    
    float total = 0;
    
    for (int i = 0; i < otus.size(); i ++)
        if (abds[i] > 0){
            
            string a_otu = Check_FMS_OTU(otus[i]);
            float cp_no = 1.0;
                                      
            if (Id_hash.count(a_otu) != 0){
                otu_abd[a_otu] = abds[i] / cp_no;
                total += otu_abd[a_otu];
                }                
            }    
    //norm
    total /= 100.0;
    int mapped_otu_count = 0;
    for (hash_map<string, float, std_string_hash>::iterator miter = otu_abd.begin(); miter != otu_abd.end(); miter ++){
                         
                          if (Id_hash.count(miter->first) != 0){
                             Abd[Id_hash[miter->first]] = miter->second;
                             mapped_otu_count ++;
                             }
                          }
    return mapped_otu_count;
    }


float _FMS_Comp_Tree::Calc_FMS_sim(float * Abd_1, float * Abd_2){
      
      float Reg_1[REG_SIZE];
      float Reg_2[REG_SIZE];
      
      float total = 0;
      
      for(int i = 0; i < OrderN; i++){
              
              int order_1 = Order_1[i];
              int order_2 = Order_2[i];
              int order_d = Order_d[i] + REG_SIZE;
              
              float dist_1 = 1- Dist_1[i];
              float dist_2 = 1- Dist_2[i];
              
              dist_1 = (dist_1 < 0) ? MIN_DIST : dist_1;
              dist_2 = (dist_2 < 0) ? MIN_DIST : dist_2;
              
              float c1_1;
              float c1_2;
              
              float c2_1;
              float c2_2;
                                
              if (order_1 >= 0){
                          
                          c1_1 = Abd_1[order_1];
                          c1_2 = Abd_2[order_1];
                          
                          }
              else {
                   c1_1 = Reg_1[order_1 + REG_SIZE];
                   c1_2 = Reg_2[order_1 + REG_SIZE];
                   }
              
              if (order_2 >= 0){
                          
                          c2_1 = Abd_1[order_2];
                          c2_2 = Abd_2[order_2];
                          
                          }
              else {
                   c2_1 = Reg_1[order_2 + REG_SIZE];
                   c2_2 = Reg_2[order_2 + REG_SIZE];
                   }
              //min
              float min_1 = (c1_1 < c1_2)?c1_1:c1_2;
              float min_2 = (c2_1 < c2_2)?c2_1:c2_2;
              
              total += min_1;
              total += min_2;
              
              //reduce
              Reg_1[order_d] = (c1_1 - min_1) * dist_1 + (c2_1 - min_2) * dist_2;
              Reg_2[order_d] = (c1_2 - min_1) * dist_1 + (c2_2 - min_2) * dist_2;
              
              }
      
        //total /= 100.0; //scale 0-1
        total = (total > 1.0) ? 1 : total;
        total = (total < 0.0) ? 0 : total;
        float abd_m=0,abd_n=0;
	int len_m=99322,len_n=99322;
	    for(int i=0;i<len_m;i++)
		    abd_m+=Abd_1[i];
	    for(int j=0;j<len_n;j++)
		    abd_n+=Abd_2[j];
            total=2*total/(abd_m+abd_n);
        return total;
      }
      
#endif
