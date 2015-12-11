%=======================================================================
%
%     getnuprox
%     function [r,d,t,ye,mue,muhat,ynu,znu,tnu,eta,...
%       yedot,edot,LN,LE] = getnuprox(filename,n);
%
%=======================================================================

      function [r,d,t,ye,mue,muhat,ynu,znu,tnu,eta,...
        yedot,edot,LN,LE] = getnuprox(filename,n);

%     function [r,d,t,ye,mue,muhat,ynu,znu,tnu,eta,...
%       yedot,edot,LN,LE] = getnuprox(filename,n);

%.....open file.........................................................
      fid = fopen(filename);

%.....skip header.......................................................
%     disp('reading nuprox')
      header = readln(fid,3);

%.....read hydro block..................................................
      block = readln(fid,n);
      [i,j] = find(block == 42.*ones(size(block)));
      if (size(i)>0)
        block(i,j) = 48;
      end

% SCW quick hack to avoid -Infinity
      [i,j] = find(block == 73.*ones(size(block)));
      if (size(i)>0)
        block(i,j) = 48;
	block(i,j+1) = 48;
	block(i,j+2) = 48;
	block(i,j+3) = 48;
	block(i,j+4) = 48;
	block(i,j+5) = 48;
	block(i,j+6) = 48;
	block(i,j+7) = 48;
      end

      block = str2num(char(block));
      r = block(:,3);
      d = block(:,4);
      t = block(:,5);
      ye = block(:,6);
      mue = block(:,7);
      muhat = block(:,8);
      ynu(:,1) = block(:,9);
      ynu(:,2) = block(:,10);
      znu(:,1) = block(:,11);
      znu(:,2) = block(:,12);
      tnu(:,1) = block(:,13);
      tnu(:,2) = block(:,14);
      eta(:,1) = block(:,15);
      eta(:,2) = block(:,16);
      yedot(:,1) = block(:,17);
      yedot(:,2) = block(:,18);
      edot(:,1) = block(:,19);
      edot(:,2) = block(:,20);
      LN(:,1) = block(:,21);
      LN(:,2) = block(:,22);
      LE(:,1) = block(:,23);
      LE(:,2) = block(:,24);

      fclose(fid);

%=======================================================================
