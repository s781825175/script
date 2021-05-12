#coding=utf-8
import pymysql
import os,sys
import re
from optparse import OptionParser
import time
usage="Annokb.py -f xxx.csv -t snv(indel) -o xxx.csv"
parser=OptionParser(usage=usage)
parser.add_option("-f","--file",dest="infile",help="input variant csv format for snp or indel")
parser.add_option("-t","--type",dest="type",help="chose variant type eg.snv,indel")
parser.add_option("-o","--output",dest="outfile",help="output file,*.csv")
(options,args)=parser.parse_args()

def Indel(cur):
	try:
		IN1=open(options.infile,'r')
		OUT=open(options.outfile,'w')
	except Exception as e:
		raise e
	title=IN1.readline().strip()
	titleline=title.replace("\"","").split(",")
	print>>OUT,title,",","KBmark",",","Catalog",",","domain",",","range",",","Length",",","Source",",","Function",",","KBmark-more"
	for line in IN1:
		line=line.strip().replace("\"","")
		arr=line.split(",")
		# if (arr[5]!="exonic" and ("deletion" not in arr[8] or "insertion" not in arr[8])) or "UNKNOWN" in arr[9]:
			# continue
		linedict=dict(zip(titleline,arr))
		# arr2=arr[9].split(",")
		# gene=arr[0].replace("\"","")
		# exon=arr[3].replace("\"","")
		# protein=arr[5].replace("\"","")
		aa_ref=[]
		site_ref=[]
		aa_var=[]
		count=0
		count2=0
		trunk=0
		if 'fs' in linedict["AA"]:
			in_ref=re.findall("p.([A-Z=*])(\d+)fs",linedict["AA"])
			if len(in_ref)>0:
				aa_ref.append(in_ref[0][0])
				site_ref.append(in_ref[0][1])
				#cur.execute("select * from all_ah_hapknow where gene=%s and variant like '%%fs%%' and variant like %s" ,('%'+linedict["gene"]+';'+'%','%'+in_ref[0][1]+'%')) #remember variant expression
				cur.execute("select * from all_ah_hapknow where gene=%s and variant like '%%fs%%' and variant like %s",(linedict["gene"],"%"+in_ref[0][1]+"%")) #remember variant expression
				all_row=cur.fetchall()
				
				for row in all_row:
					aa_row=re.findall("([A-Z=*])(\d+)([A-Z=*])fs",row[10],re.I)
					if len(aa_row)<1 or len(aa_row[0])<3:
						aa_row=re.findall("([A-Z=*])(\d+)fs",row[9],re.I)
					if len(aa_row)<1 or len(aa_row[0])<2:
						continue
					for i in range(len(aa_row)):
						if aa_row[i][0]==aa_ref[0] and aa_row[i][1]==site_ref[0]:
							count=1
						if (aa_row[i][0]=="X" and aa_row[i][1]==site_ref[0] and aa_ref[0]=="*") or (aa_row[i][0]=="*" and aa_row[i][1]==site_ref[0] and aa_ref[0]=="X"):
							count2=1
		if 'delins' in linedict["AA"]:
			in_ref=re.findall("p.([A-Z=*])(\d+)delins([A-Z=*]+)",linedict["AA"])
			if len(in_ref)<1:
				in_ref=re.findall("p.([A-Z=*])(\d+)delins([A-Z=*])",linedict["AA"])
			if len(in_ref)>0:
				cur.execute("select * from all_ah_hapknow where gene=%s and variant like '%%delins%%' and variant like %s ",(linedict["gene"],"%"+in_ref[0][1]+"%"))
				all_row=cur.fetchall()
				aa_ref.append(in_ref[0][0])
				site_ref.append(in_ref[0][1])
				aa_var.append(in_ref[0][2])				
				for row in all_row:
					#print row[9]
					#row[9]=row[9].replace("p.","")
					aa_row=re.findall("([A-Z=*])(\d+)delins([A-Z=*]+)",row[10],re.I)
					if len(aa_row)<1 or len(aa_row[0])<3:
							continue
					for i in range(len(aa_row)):
						if aa_row[i][0]==aa_ref[0] and aa_row[i][1]==site_ref[0] and aa_row[i][2]==aa_var[0]:
						#ll=Print(row)
							count=1
						if ((aa_row[i][0] in "X*" and aa_ref[0] in "X*") and aa_row[i][1]==site_ref[0] and aa_row[i][2]==aa_var[0]) or (aa_row[i][0]==aa_ref[0] and aa_row[i][1]==site_ref[0] and (aa_row[i][2] in "X*" and aa_var[0] in "X*")):
							trunk=1
						if aa_row[i][0]==aa_ref[0] and aa_row[i][1]==site_ref[0] and (aa_row[i][2]!=aa_var[0] and trunk==0):
							count2=1
						if trunk==1 and count==0:
							count=1
			else:
				in_ref=re.findall("p.([A-Z=*])(\d+)_([A-Z=*])(\d+)delins([A-Z=*]+)",linedict["AA"])
				if len(in_ref)<1:
					in_ref=re.findall("p.([A-Z=*])(\d+)_([A-Z=*])(\d+)delins([A-Z=*])",linedict["AA"])
				if len(in_ref)>0:
					cur.execute("select * from all_ah_hapknow where gene=%s and variant like '%%delins%%' and variant like %s ",(linedict["gene"],"%"+in_ref[0][1]+"%"))
					all_row=cur.fetchall()
					aa_ref.append(in_ref[0][0])
					aa_ref.append(in_ref[0][2])
					site_ref.append(in_ref[0][1])
					site_ref.append(in_ref[0][3])
					aa_var=in_ref[0][4]
					for row in all_row:
						#print row[9]
						#row[9]=row[9].replace("p.","")
						aa_row=re.findall("([A-Z=*])(\d+)_([A-Z=*])(\d+)delins([A-Z=*]+)",row[10],re.I)
						if len(aa_row)<1 or len(aa_row[0])<5:
							continue
						for i in range(len(aa_row)):
							if aa_row[i][0]==aa_ref[0] and aa_row[i][1]==site_ref[0] and aa_row[i][2]==aa_ref[1] and aa_row[i][3]==site_ref[1] and aa_row[0][4]==aa_var[1]:
							#ll=Print(row)
								count=1
							if ((aa_row[i][0] in "X*" and aa_ref[0] in "X*") and aa_row[i][1]==site_ref[0] and aa_row[i][2]==aa_ref[1] and aa_row[i][3]==site_ref[1] and aa_row[0][4]==aa_var[1]) or (aa_row[i][0]==aa_ref[0] and aa_row[i][1]==site_ref[0] and (aa_row[i][2] in "X*" and a_ref[1] in "X*") and aa_row[i][3]==site_ref[1] and aa_row[0][4]==aa_var[1]):
			#print(in_ref)		
								trunk=1
							if aa_row[i][0]==aa_ref[0] and aa_row[i][1]==site_ref[0] and aa_row[i][2]==aa_ref[1] and aa_row[i][3]==site_ref[1] and (aa_row[0][4]!=aa_var[1] and trunk==0):
								count2=1
							if trunk==1 and count==0:
								count=1
		if 'ins' not in linedict["AA"] and 'del' in linedict["AA"]:
			
			in_ref=re.findall("p.(\d+)_(\d+)del",linedict["AA"])
			if len(in_ref)<1:
				in_ref=re.findall("p.[A-Z=*](\d+)_[A-Z=*](\d+)del",linedict["AA"])
			if(len(in_ref)>0):
				cur.execute("select * from all_ah_hapknow where gene=%s and variant not like '%%ins%%' and variant like '%%del%%' and variant like %s " ,(linedict["gene"],"%"+in_ref[0][1]+"%"))
				all_row=cur.fetchall()
				site_ref.append(in_ref[0][0])
				site_ref.append(in_ref[0][1])
				for row in all_row:
					aa_row=re.findall("([A-Z=*])(\d+)_([A-Z=*])(\d+)del",row[10],re.I)
					if len(aa_row)>0 and len(aa_row[0])==4:
						for i in range(len(aa_row)):
							if aa_row[i][1]==site_ref[0] and aa_row[i][3]==site_ref[1]:
							#ll=Print(row)
								count=1
							# if ((aa_row[i][0] in "X*" and site_ref[1] in "X*") and aa_row[i][1]==site_ref[1]):
								# count=1
					else:
						aa_row=re.findall("([A-Z=*])(\d+)del",row[10],re.I)
						if len(aa_row)<1 or len(aa_row[0])<2:
							continue
						for i in range(len(aa_row)):
							if aa_row[i][0]==site_ref[0] and aa_row[i][1]==site_ref[1]:
							#ll=Print(row)
								count=1
							if (aa_row[i][0] in "X*" and site_ref[0] in "X*") and aa_row[i][1]==site_ref[1]:
								count=1
						else:
							aa_row=re.findall("(\d+)_(\d+)del",row[10],re.I)
							if len(aa_row)<1 or len(aa_row[0])<2:
								continue
							for i in range(len(aa_row)):
								if aa_row[i][0]==site_ref[0] and aa_row[i][1]==site_ref[1]:
									count=1
							
		print>>OUT,line,",",count,","," ",","," ",","," ",","," ",","," ",","," ",",",count2
		OUT.flush()	
	IN1.close()
	OUT.close()
def Snp(cur):
	try:
		IN1=open(options.infile,'r')
		OUT=open(options.outfile,'w')
	except Exception as e:
		raise e
	title=IN1.readline().strip()
	titleline=title.replace("\"","").split(",")
	print>>OUT,title,",","KBmark",",","Catalog",",","domain",",","range",",","Length",",","Source",",","Function",",","KBmark-more"
	for line in IN1:
		
		line=line.strip().replace("\"","")
		arr=line.split(",")
		linedict=dict(zip(titleline,arr))
		# if ((arr[5]!="exonic" and "SNV" not in arr[8]) or "UNKNOWN" in arr[9]):
			# continue
		# gene=arr[0].replace("\"","")
		# exon=arr[3].replace("\"","")
		# protein=arr[5].replace("\"","")
		aa_ref=[]
		site_ref=[]
		aa_var=[]
		#cur.execute("select * from all_ah_hapknow where Gene like %s and Variant not like '%%ins%%' and Variant not like '%%del%%' and Variant not like '%%\_%%' ",(linedict["gene"]+';'))
		#cur.execute("select * from db_collects where gene=%s and site=%s and variant=%s",(gene,exon,protein))
		in_ref=re.findall("p.([A-Z=*])(\d+)([A-Z=*])",linedict["AA"],re.I)
		count=0
		count2=0
		trunk=0
		if len(in_ref)>0:
			gg=str(in_ref[0][1])
			
			cur.execute("select * from all_ah_hapknow where gene=%s and variant like %s and variant like %s",(linedict["gene"],"%"+gg+"%",'%'+in_ref[0][1]+'%'))
			#cur.execute("select * from all_ah_hapknow where Gene like %s ",('%'+linedict["gene"]+';'+'%'))
			#print "select * from all_ah_hapknow where Gene like %s and Variant like %s and Variant like %s"%('%'+linedict["gene"]+';'+'%',"%"+gg+"%",'%'+in_ref[0][1]+'%')
			
			aa_ref.append(in_ref[0][0])
			site_ref.append(in_ref[0][1])
			aa_var.append(in_ref[0][2])
			all_row=cur.fetchall()
			for row in all_row:
				aa_row=re.findall("([A-Z*=])(\d+)([A-Z*=])",row[10],re.I)
				if len(aa_row)<1 or len(aa_row[0])<3:
					continue
				for i in range(len(aa_row)):
					if aa_row[i][0]==aa_ref[0] and aa_row[i][1]==site_ref[0] and (aa_row[i][2]==aa_var[0] or (aa_row[i][2]=="=" and aa_var[0]==aa_ref[0])):
						count=1
					if ((aa_row[i][0]=="X" and aa_ref[0]=="*") or (aa_row[i][0]=="*" and aa_ref[0]=="X")) and aa_row[i][1]==site_ref[0] and (aa_row[i][2]==aa_var[0] or (aa_row[i][2]=="=" and aa_var[0]==aa_ref[0])):
						trunk=1
					if aa_row[i][1]==site_ref[0] and aa_row[i][0]==aa_ref[0] and ((aa_row[i][2]=="*" and aa_var[0]=="X") or (aa_row[i][2]=="X" and aa_var[0]=="*")):
						trunk=1
					if aa_row[i][0]==aa_ref[0] and aa_row[i][1]==site_ref[0] and (aa_row[i][2]!=aa_var[0] and trunk==0):
						count2=1
					if ((aa_row[i][0]=="X" and aa_ref[0]=="*") or (aa_row[i][0]=="*" and aa_ref[0]=="X")) and aa_row[i][1]==site_ref[0] and (aa_row[i][2]!=aa_var[0] and trunk==0):
						count2=1
					if trunk==1 and count==0:
						count=1
					
			print>>OUT,line,",",count,","," ",","," ",","," ",","," ",","," ",","," ",",",count2
		else:
			count=0
			print>>OUT,line,",",count,","," ",","," ",","," ",","," ",","," ",","," ",",",count2
	IN1.close()
	OUT.close()
# def Print(row_list):
	# row_ll=list(row_list)
	# row_ll[0]=str(row_ll[0])
	# ll="\t".join(row_ll)
	# return ll

def main():
	conn=pymysql.connect(host='192.168.1.10',port=3306,user='hapknow',passwd='Haplox2017!',db='hapknow')
	cursor=conn.cursor()
	print time.asctime( time.localtime(time.time()) )
	if options.type=="indel":
		Indel(cursor)
	if options.type=="snv":
		Snp(cursor)
	print time.asctime( time.localtime(time.time()) )
#	if options.type=="fusion":
#		Fusion(cursor)
#	if options.type=="cnv":
#		Cnv(cursor)


if __name__ == "__main__":		
	main()
