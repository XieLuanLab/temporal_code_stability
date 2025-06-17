function pressZeroEnter(HObj, event)
  commandwindow
  import java.awt.*;
  import java.awt.event.*;
  rob = Robot;
  pause(0.01)
  rob.keyPress(KeyEvent.VK_NUMPAD0)
  rob.keyRelease(KeyEvent.VK_NUMPAD0)
  pause(0.01)
  rob.keyPress(KeyEvent.VK_ENTER)
  rob.keyRelease(KeyEvent.VK_ENTER)
  pause(0.01)
end