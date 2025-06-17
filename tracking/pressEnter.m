function pressEnter(HObj, event)
  import java.awt.*;
  import java.awt.event.*;
  rob = Robot;
  pause(0.01)
  rob.keyPress(KeyEvent.VK_ENTER)
  rob.keyRelease(KeyEvent.VK_ENTER)
  pause(0.01)
end