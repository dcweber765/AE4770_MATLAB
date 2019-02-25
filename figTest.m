figure(1) = openfig('liftDragTrade.fig');
   figure(2) = openfig('rangeTrade.fig');
   figure(3) = openfig('edurTrade.fig');

   L3 = findobj(2,'type','line');
   L4 = findobj(3,'type','line');
  copyobj(L3,findobj(1,'type','axes'));
  copyobj(L4,findobj(1,'type','axes'));