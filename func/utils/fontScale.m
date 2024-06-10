function fontScale(scale)
% function fontScale(scale)

H = gcf;
allText   = findall(H, 'type', 'text');
allAxes   = findall(H, 'type', 'axes');
allFont   = [allText; allAxes];
fontSize = get(allFont,'FontSize');

newFontSize = LocalScale(fontSize, scale, 2);
set(allFont,{'FontSize'},newFontSize);

function newArray = LocalScale(inArray, scale, minValue)
n = length(inArray);
newArray = cell(n,1);
for k=1:n
  newArray{k} = max(minValue,scale*inArray{k}(1));
end