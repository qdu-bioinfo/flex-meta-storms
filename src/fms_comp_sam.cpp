// Updated at Sept 30, 2022
// Updated by Mingqian Zhang
// Bioinformatics Group, College of Computer Science & Technology, Qingdao University
// Code by: Mingqian Zhang, Xiaoquan Su
// Flex Meta-Storms version 1.0

#include <iostream>
#include <fstream>
#include <sstream>
#include<string>
#include<vector>
#include<algorithm>
#include<unordered_map>

#include <sys/stat.h>
#include <unistd.h>
#include <omp.h>
#include <string.h>
#include <stdlib.h>

#include "fms_func.h"
#include "fms_comp.h"

using namespace std;

char Ref_db;

string Listfilename;
string Listprefix;
string Queryfile1;
string Queryfile2;

string Tablefilename;
string Outfilename;
string Targetfilename;
string Markerfilename;

int Coren = 0;

bool Is_cp_correct; 
bool Is_sim; //true: sim, false: dist;
bool Is_heatmap;
int Cluster = 2;

int Mode = 1; //default is FMS;

 
int printhelp(){
    
    cout << "Version : " << Version << endl;
    cout << "\tCompute the Flex Meta-Storms similarity/distance among samples" << endl;
    cout << "Usage: " << endl;
    cout << "FMS-comp-taxa [Option] Value" << endl;
    cout << "Options: " << endl;
    
    cout << "\t[Input options, required]" << endl;
    cout << "\t  -T (upper) Input OTU count table (*.OTU.Abd) for multi-sample comparison" << endl;
	cout << "\t  -M Enter biomarker path (.xls)" << endl;

    cout << "\t  -K Calculate Exact marker distance" << endl;
	cout << "\t  or " << endl;
    cout << "\t  -L Calculate Flex Meta-Storms distance" << endl;
    
    cout << "\t[Output options]" << endl;
    cout << "\t  -o Output file, default is to output on screen" << endl;
    cout << "\t  -d Output format, distance (T) or similarity (F), default is T" << endl;
    cout << "\t[Other options]" << endl;
    cout << "\t  -h Help" << endl;
    
    exit(0);
    
    return 0;
    
    };
    
int Parse_Para(int argc, char * argv[]){
    
    Ref_db = 'G';
    
    Coren = 0;
    Mode = 1; //default is FMS;
    
    Is_cp_correct = true;
    Is_sim = false;
    Is_heatmap = false;
    
    int i = 1;
      
      if (argc ==1) 
		printhelp();
      
      while(i<argc){
         if (argv[i][0] != '-') {
                           cerr << "Argument # " << i << " Error : Arguments must start with -" << endl;
                           exit(0);
                           };           
         switch(argv[i][1]){
                            case 'D': Ref_db = argv[i+1][0]; break;
                            case 'T': Tablefilename = argv[i+1]; break;
							case 'K': Targetfilename = ".exact_marker"; Mode = 0; i -= 1; break;
							case 'L': Targetfilename = ".target_marker"; Mode = 1; i -= 1; break;
                            case 'o': Outfilename = argv[i+1]; break;
							case 'M': Markerfilename = argv[i + 1]; break;
			    case 'd': if ((argv[i + 1][0] == 'f') || (argv[i + 1][0] == 'F')) Is_sim = true; break;
                            case 't': Coren = atoi(argv[i+1]); break;         
                            case 'h': printhelp(); break;
                            default : cerr << "Error: Unrec argument " << argv[i] << endl; printhelp(); break; 
                            }
         i+=2;
         }
    
    int max_core_number = sysconf(_SC_NPROCESSORS_CONF);
    
    if ((Coren <= 0) || (Coren > max_core_number)){
                    
                    Coren = max_core_number;
                    } 
    if (Cluster <= 0){
                cerr << "Warning: cluster number must be larger than 0, change to default (2)" << endl;
                }
    }

vector<string> bio_marker;
void Read_bio_marker(string path1)
{
	ifstream ifs_bio_marker;
	ifs_bio_marker.open(path1, ios::in);
	if (!ifs_bio_marker)
	{
		cerr << "Error: Cannot open file : " << path1 << endl;
		return;
	}
	string temp;
	vector<string> line;
	getline(ifs_bio_marker, temp);
	while (getline(ifs_bio_marker, temp))
	{
		line = split(temp, "\t");
		bio_marker.emplace_back(line[0]);
		line.clear();
	}
		
}

unordered_map<string, vector<string>> similarity_table;
void save_similarity_table()
{
	string PATH;
	PATH = Check_FMS_Env();
	PATH += "/databases/gg_13/gg_13_approximate_members.tab";
	ifstream ifs2;
	ifs2.open(PATH, ios::in);
	vector<string> sim_line;
	vector<string> sim_temp;
	string temp;
	while (getline(ifs2, temp))
	{
		sim_line = split(temp, "\t");
		temp = sim_line[0];
		if (sim_line.size() > 1)
		{
			for (int i = 1; i < sim_line.size(); i++)
				sim_temp.emplace_back(sim_line[i]);
		}
		else
			sim_temp.emplace_back("00");
		similarity_table.insert(pair < string, vector<string>>(temp, sim_temp));
		sim_temp.clear();
	}
	ifs2.close();
}


vector<string> sample_otu_sum;
void save_sample_otu_sum(string Tablefilename)
{
	ifstream ifs_sample_otu;
	ifs_sample_otu.open(Tablefilename, ios::in);
	string tempt;
	getline(ifs_sample_otu, tempt);
	sample_otu_sum = split(tempt, "\t");
	ifs_sample_otu.close();
}


unordered_map<string, float> map_local_otu;
unordered_map<string, float> map_key_otu;
void exact_marker(string Tablefilename)
{
	ifstream ifs_talbe_file;
	ifs_talbe_file.open(Tablefilename, ios::in);
	ofstream ofs_exact_marker;
	Targetfilename = Outfilename + Targetfilename;
	ofs_exact_marker.open(Targetfilename, ios::out);
	vector<vector<string>> exact_marker;
	vector<float> single_sample;
	vector<string> vec_temp;
	string temp;
	getline(ifs_talbe_file, temp);
	int times = 0;
	cout << endl;
	cout << "Calculate exact markers distance :" << endl;
	cout << "Number of exact markers is : " << bio_marker.size() << endl;
	while (getline(ifs_talbe_file, temp))
	{
		vec_temp = split(temp, "\t");
		for (int i = 0; i < vec_temp.size(); i++)
		{
			single_sample.push_back(stringtonum<float>(vec_temp[i]));
		}
		for (int i = 0; i < sample_otu_sum.size(); i++)
		{
			for (int j = 0; j < bio_marker.size(); j++)
			{
				if (sample_otu_sum[i] == bio_marker[j])
				{
					map_key_otu[bio_marker[j]] = single_sample[i];
				}
			}
		}
		if (times == 0)
		{
			ofs_exact_marker << "sample_id" << "\t";
			for (auto pointer = map_key_otu.begin(); pointer != map_key_otu.end(); pointer++)
			{
				if (times == 0)
					ofs_exact_marker << pointer->first << "\t";
			}
			ofs_exact_marker << endl;
		}
		ofs_exact_marker << vec_temp[0] << "\t";
		for (auto pointer = map_key_otu.begin(); pointer != map_key_otu.end(); pointer++)
		{
			ofs_exact_marker << pointer->second << "\t";
		}
		ofs_exact_marker << endl;

		for (auto pointer = map_key_otu.begin(); pointer != map_key_otu.end(); pointer++)
		{
			pointer->second = 0;
		}
		single_sample.clear();
		vec_temp.clear();
		times++;
	}
}
void target_marker(string Tablefilename)
{
	ifstream ifs_talbe_file;
	ifs_talbe_file.open(Tablefilename, ios::in);
	ofstream  ofs_target;
	Targetfilename = Outfilename + Targetfilename;
	ofs_target.open(Targetfilename, ios::out);
	vector<vector<string>> exact_marker;
	vector<float> single_sample;
	vector<string> vec_temp;
	string temp;
	getline(ifs_talbe_file, temp);
	int times = 0;
	while (getline(ifs_talbe_file, temp))
	{
		vec_temp = split(temp, "\t");
		for (int i = 0; i < vec_temp.size(); i++)
		{
			single_sample.push_back(stringtonum<float>(vec_temp[i]));
		}
		for (int i = 0; i < sample_otu_sum.size(); i++)
		{
			for (int j = 0; j < bio_marker.size(); j++)
			{
				if (sample_otu_sum[i] == bio_marker[j])
				{
					map_local_otu[bio_marker[j]] = single_sample[i];
				}
			}
		}

		for (int i = 0; i < bio_marker.size(); i++)
		{
			for (int j = 1; j < sample_otu_sum.size(); j++)
			{
				for (int m = 0; m < similarity_table[bio_marker[i]].size(); m += 2)
				{
					if (sample_otu_sum[j] == similarity_table[bio_marker[i]][m])
					{
						vector<string> ::iterator t;
						t = find(bio_marker.begin(), bio_marker.end(), sample_otu_sum[j]);
						if (t == bio_marker.end()) {
							if (map_local_otu.count(sample_otu_sum[j]) == 0)
								map_local_otu.insert(pair<string, float>(sample_otu_sum[j], stringtonum<float>(similarity_table[bio_marker[i]][m + 1]) * single_sample[j]));
							else
							{
								if (map_local_otu[sample_otu_sum[j]] < stringtonum<float>(similarity_table[bio_marker[i]][m + 1]) * single_sample[j])
								{
									map_local_otu[sample_otu_sum[j]] = stringtonum<float>(similarity_table[bio_marker[i]][m + 1]) * single_sample[j];
								}
							}
						}


					}
				}
			}
		}
		if (times == 0)
		{
			ofs_target << "sample_id" << "\t";
			for (auto pointer = map_local_otu.begin(); pointer != map_local_otu.end(); pointer++)
			{
				ofs_target << pointer->first << "\t";
			}
			ofs_target << endl;
			cout << endl;
			cout << "Calculate flex meta-storms distance :" << endl;
			cout << "Number of exact markers is : " << bio_marker.size() << endl;
			cout << "Number of approximate markers is : " << map_local_otu.size()-bio_marker.size() << endl;
			cout << "Number of target markers is : " << map_local_otu.size() << endl;
		}
		ofs_target << vec_temp[0] << "\t";
		for (auto pointer = map_local_otu.begin(); pointer != map_local_otu.end(); pointer++)
		{
			ofs_target << pointer->second << "\t";

		}
		for (auto pointer = map_local_otu.begin(); pointer != map_local_otu.end(); pointer++)
		{
			pointer->second = 0;
		}
		ofs_target << endl;
		single_sample.clear();
		vec_temp.clear();
		times++;
	}
}

void Output_Matrix(const char * outfilename, int n, vector <float> * sim_matrix, bool is_sim, vector <string> sam_name){
	
     ofstream outfile(outfilename, ofstream::out);
     if (!outfile){
                   cerr << "Error: Cannot open output file : " << outfilename << endl;
                   return; 
                   }
     
     //Label
     for(int i = 0; i < n; i ++)
             outfile << "\t" << sam_name[i];
     outfile << endl;
     
     for(int i = 0; i < n; i ++){
             outfile << sam_name[i];
             
             for (int j = 0; j < n; j ++){                
                 
                 long ii = (i <= j) ? i : j;
                 long jj = (i >= j) ? i : j;
                 long p = ii * (long) n + jj - (1 + ii + 1) * (ii + 1) / 2;
                                  
                 if (is_sim){                 
                    if (ii == jj) outfile << "\t" << 1.0;
                    else outfile << "\t" << (*sim_matrix)[p];
                    }
                 else {
                      if (ii == jj) outfile << "\t" << 0.0;
                      else outfile << "\t" << 1.0 - (*sim_matrix)[p];
                      }
                 }
             outfile << endl;
             }
                 
     outfile.close();
     outfile.clear();
     }


void Multi_Comp_Table(_Table_Format abd_table){
	
    _FMS_Comp_Tree comp_tree(Ref_db);
    int file_count = abd_table.Get_Sample_Size();
         
    //load abd
    float **Abd = new float * [file_count];
    for (int i = 0; i < file_count; i ++){
        Abd[i] = new float [comp_tree.Get_LeafN()];
        cout << comp_tree.Load_FMS_abd(&abd_table, Abd[i], i, Is_cp_correct) << " OTUs in file " << i + 1 << endl;
        }
    
    cout << file_count << " files loaded" << endl;
    
    //make order
    vector <int> order_m;
    vector <int> order_n;
    long iter = 0;
    for (int i = 0; i < file_count - 1; i ++)
        for (int j = i + 1; j < file_count; j ++){            
            order_m.push_back(i);
            order_n.push_back(j);
            iter ++;
            }
        
    vector <float>  sim_matrix;
    for (long i = 0; i < iter; i ++)
        sim_matrix.push_back(0);
    
    //openmp
        
    omp_set_num_threads(Coren);
    
    #pragma omp parallel for schedule(dynamic, 1)
    for (long i = 0; i < iter; i ++){
        
        long m = order_m[i];
        long n = order_n[i];
        long p = m * (long) file_count + n - (1 + m + 1) * (m + 1) / 2;
        
        sim_matrix[p] = comp_tree.Calc_FMS_sim(Abd[m], Abd[n]);
        }
    
    Output_Matrix(Outfilename.c_str(), file_count, &sim_matrix, Is_sim, abd_table.Get_Sample_Names());
    
    for (int i = 0; i < file_count; i ++)
        delete [] Abd[i];
               
     };

int main(int argc, char * argv[]){

    Parse_Para(argc, argv);     
	Read_bio_marker(Markerfilename);
	save_similarity_table();
	save_sample_otu_sum(Tablefilename);
	
    switch (Mode){
           case 0: {
			   exact_marker(Tablefilename);
			   
			   _Table_Format table(Targetfilename.c_str());
			   Multi_Comp_Table(table);
			   break;
		   }
           case 1:{
			   target_marker(Tablefilename);

               _Table_Format table(Targetfilename.c_str());
               Multi_Comp_Table(table); 
               break;
                   }
           default: break;
           }
    return 0;
    }
