path_to_matrix = './'
path_to_cpgPos = './'
name = 'test_chr1_1287967_1288117';

data1 = load(strcat(path_to_matrix,'/',name,'_me.txt'));
data2 = load(strcat(path_to_matrix,'/',name,'_unme.txt'));
merge = data1-data2;

merge_v = merge(:);
colorz = zeros(length(merge_v),3);
colorz(merge_v == -1,:) = repmat([0.68,0.92,1],sum(merge_v == -1),1); % unmethylated CpGs are labeled as light blue
colorz(merge_v == 0,:) = repmat([0.8,0.8,0.8],sum(merge_v == 0),1); % CpGs without signal are labeled as grey
colorz(merge_v == 1,:) = repmat([0.6,0.2,0],sum(merge_v == 1),1); % methylated CpGs are labeled as dark red


x = [];
for i = 1:size(merge,2)
    x = [x,i*ones(size(merge,1),1)];
end

tmp = 1:size(merge,1);
y = repmat(tmp,1,size(merge,2));


figure(1);
plot([0.5,size(data1,2)+0.5],ones(1,2));
for i = 2:size(data1,1)
    hold on;plot([0.5,size(data1,2)+0.5],repmat(i,1,2),'k');
end
hold on; scatter(x(:),y(:),30,colorz,'o','filled','markeredgecolor','black');
xlim([0,size(data1,2)+1]);
ylim([0,size(data1,1)+1]);title(strrep(name, '_', '-'));
hold off;


%% with distance between CpGs
out_split = strsplit(name,'_');
chrom = out_split{end-2};
start = str2num(out_split{end-1});
en = str2num(out_split{end});
cpgpos = readtable(strcat(path_to_cpgPos,'/cpgpos_',chrom,'.pos'),...
    'ReadVariableNames',false,'filetype','text');

posw = cpgpos{cpgpos{:,2}>= start & cpgpos{:,2}< en,2};
relative_posw = posw - start;

x2 = [];
for i = 1:size(merge,2)
    x2 = [x2,relative_posw(i)*ones(size(merge,1),1)];
end

figure(2);
plot([-5,en-start+4],ones(1,2),'k');
for i = 2:size(data1,1)
    hold on;plot([-5,en-start+4],repmat(i,1,2),'k');
end
hold on; scatter(x2(:),y(:),30,colorz,'o','filled','markeredgecolor','black');
xlim([-8,en-start+7]);
ylim([0,size(data1,1)+1]);title(strrep(name, '_', '-'));
hold off;