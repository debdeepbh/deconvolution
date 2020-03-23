function demo_save(fignum, filename, sizecomma)
set(fignum,'PaperPositionMode','auto');  
pathname = 'pic/';
fullpath = strcat(pathname,filename,'.pdf');
sizeval = strcat('-S',sizecomma);	% e.g. sizecommaa = '800,200'
print(fignum,'-dpdf',sizeval,fullpath)
