%=======================================================================
%
%     getagile
%
%=======================================================================

      function [nq,ny,t,label,y,dx,s,cs,status] = getagile(filename);

%.....read header.......................................................
      fid = fopen(filename);
      [name,count] = fscanf(fid,' %*[=] file: %s');
      [n,count] = fscanf(fid,' size: %d %d',[2 1]);
      nq = n(1);
      ny = n(2)-1;
      [t,count] = fscanf(fid,' time: %g',1);
      comment = readln(fid,8);
      for iy=1:ny
        i = 16 + 15*(iy-1);
        label(iy,:) = [comment(5,i:i+14) comment(6,i:i+14)];
      end
      label = char(label);

%.....read data.........................................................
      block = readln(fid,nq);
      [i,j] = find(block == 42.*ones(size(block)));
      if (size(i)>0)
        block(i,j) = 48;
      end
      block = str2num(char(block));

%.....write into variable y.............................................
      for iy=1:ny
        y(:,iy) = block(:,1+iy);
      end

%.....skip separator....................................................
      comment = readln(fid,3);

%.....read data.........................................................
      block = readln(fid,nq);
      [i,j] = find(block == 42.*ones(size(block)));
      if (size(i)>0)
        block(i,j) = 48;
      end
      block = str2num(char(block));
      if size(block,1)==0
        dx = 0;
        s = 0;
        cs = 0;
        status = 0;
        return
      end

%.....write into variable dx............................................
      dx(:,1) = block(:,2);
      dx(:,2) = block(:,3);
      if size(block,2)>=4
        s = block(:,4);
        cs = block(:,5);
        status = block(:,6);
      else
        s = 0;
        cs = 0;
        status = 0;
      end

%.....reconstruct precision of y(:,1) from dx...........................
      x(1) = dx(1,1);
      for iq=2:nq
        x(iq) = x(iq-1) + dx(iq,1);
      end
      y(:,1) = x' + dx(:,2);
      
%=======================================================================
